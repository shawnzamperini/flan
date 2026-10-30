/**
* @file read_s3xe.cpp
*
* @brief Reads a SOLEDGE3X (S3X) simulation into a Background.
*
* SOLEDGE3X splits its grid into zones to handle X-point geometries. Each
* zone is a logically rectangular block of cells in (psi, theta, phi). In
* Flan, every SOLEDGE3X zone becomes a zone with index-space coordinates:
* x is the psi cell index (cell edges 0, 1, ..., Npsi), y the theta cell
* index (edges 0, 1, ..., Ntheta) and z the toroidal angle phi. The
* neighbors table in the mesh file gives what is across each face of each
* zone.
*
* The mesh file holds the cell corners in (R, Z), normalized as
* R = R0 + a0 * Rcorners and Z = a0 * Zcorners, and B normalized by B0,
* for a single toroidal plane (the mesh is axisymmetric). The geometry
* (positions, basis vectors and Jacobian) is computed from the corners for
* every toroidal plane of the plasma.
*/

#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "background.h"
#include "errors.h"
#include "flan_types.h"
#include "hdf5_utils.h"
#include "read_s3xe.h"
#include "vectors.h"
#include "zones.h"

namespace S3X
{
	namespace
	{
		using Vec3 = std::array<double, 3>;

		// Name of a zone's group in the SOLEDGE3X files, numbered from 1
		std::string zone_group(const int zone)
		{
			return "/zone" + std::to_string(zone + 1);
		}

		Vec3 cross(const Vec3& a, const Vec3& b)
		{
			return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
				a[0] * b[1] - a[1] * b[0]};
		}

		double dot(const Vec3& a, const Vec3& b)
		{
			return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
		}

		// a * u + b * v
		Vec3 combine(const double a, const Vec3& u, const double b,
			const Vec3& v)
		{
			return {a * u[0] + b * v[0], a * u[1] + b * v[1],
				a * u[2] + b * v[2]};
		}

		// Cylindrical unit vectors at toroidal angle phi in the Cartesian
		// frame (see phi_sign)
		Vec3 R_hat(const double phi)
			{return {std::cos(phi), phi_sign * std::sin(phi), 0.0};}

		Vec3 phi_hat(const double phi)
			{return {-std::sin(phi), phi_sign * std::cos(phi), 0.0};}

		constexpr Vec3 Z_hat {0.0, 0.0, 1.0};

		// Empty ZonedVector shaped like the zones (with their ghost cells)
		Vectors::ZonedVector<BkgFPType> zoned(const Zones::Zones& zones,
			const int dim1)
		{
			return Vectors::ZonedVector<BkgFPType>(dim1, zones.get_nx_all(),
				zones.get_ny_all(), zones.get_nz_all(), zones.get_ghost());
		}

		// Stops the run if a dataset doesn't have the expected shape
		void check_shape(const std::vector<hsize_t>& shape,
			const std::vector<hsize_t>& want, const std::string& name,
			const std::string& path)
		{
			if (shape == want) return;
			auto str = [](const std::vector<hsize_t>& s)
			{
				std::string out {"("};
				for (std::size_t n {}; n < s.size(); ++n)
					out += (n ? ", " : "") + std::to_string(s[n]);
				return out + ")";
			};
			Errors::fatal_error("Dataset " + name + " in " + path
				+ " has shape " + str(shape) + ", but should have shape "
				+ str(want) + ".");
		}

		// Reads a per-zone dataset of doubles and checks its shape
		std::vector<double> read_checked(const Hdf5::File& file,
			const std::string& name, const std::vector<hsize_t>& want)
		{
			auto ds {Hdf5::read_dataset<double>(file, name)};
			check_shape(ds.shape, want, name, file.path());
			return std::move(ds.data);
		}

		// Grid with n cells between lo and hi, returned as (centers, edges)
		std::pair<std::vector<double>, std::vector<double>> uniform_grid(
			const int n, const double lo, const double hi)
		{
			std::vector<double> centers(n), edges(n + 1);
			for (int i {}; i <= n; ++i) edges[i] = lo + i * (hi - lo) / n;
			for (int i {}; i < n; ++i)
				centers[i] = 0.5 * (edges[i] + edges[i + 1]);
			return {centers, edges};
		}
	}

	/**
	* Prints the assumption made about SOLEDGE3X that still needs to be
	* checked: which way the toroidal angle goes.
	*/
	void print_assumptions()
	{
		std::cout << "Warning! SOLEDGE3X's (R, Z, phi) is assumed to be "
			"right-handed, so phi increases clockwise seen from above. This "
			"decides the direction of the toroidal field and of the ExB "
			"drift, and still needs to be confirmed.\n";
	}

	/**
	* Reads the reference values SOLEDGE3X output is normalized by.
	*
	* Every value must be positive, otherwise the run is stopped.
	*
	* @param ref_path Path to the reference file, e.g., refParam_raptorX.h5.
	* @return The reference values.
	*/
	Reference read_reference(const std::string& ref_path)
	{
		const Hdf5::File file(ref_path);
		Reference ref {};
		const std::array<std::pair<const char*, double*>, 8> values {{
			{"/n0", &ref.n0}, {"/T0", &ref.T0}, {"/c0", &ref.c0},
			{"/tau0", &ref.tau0}, {"/phi0", &ref.phi0}, {"/R0", &ref.R0},
			{"/a0", &ref.a0}, {"/B0", &ref.B0}}};
		for (const auto& [name, value] : values)
		{
			*value = Hdf5::read_scalar<double>(file, name);
			if (!(*value > 0.0))
			{
				Errors::fatal_error(std::string(name) + " in " + ref_path
					+ " is " + std::to_string(*value) + ", but should be "
					"positive.");
			}
		}
		return ref;
	}

	/**
	* Prints the reference values used to unnormalize SOLEDGE3X output.
	*
	* @param ref The reference values.
	*/
	void print_reference(const Reference& ref)
	{
		std::cout << "SOLEDGE3X reference values: n0 = " << ref.n0
			<< " m-3, T0 = " << ref.T0 << " eV, c0 = " << ref.c0
			<< " m/s, tau0 = " << ref.tau0 << " s, phi0 = " << ref.phi0
			<< " V\n";
	}

	/**
	* Checks that the reference file and the mesh file agree on R0, a0 and
	* B0, which both have, and stops the run if not. A mismatch means they
	* are from different SOLEDGE3X cases.
	*
	* @param ref The reference values.
	* @param mesh The mesh.
	*/
	void check_reference(const Reference& ref, const Mesh& mesh)
	{
		const std::array<std::tuple<const char*, double, double>, 3> pairs {{
			{"R0", ref.R0, mesh.R0}, {"a0", ref.a0, mesh.a0},
			{"B0", ref.B0, mesh.B0}}};
		for (const auto& [name, a, b] : pairs)
		{
			if (std::abs(a - b) > 1e-6 * std::abs(b))
			{
				Errors::fatal_error(std::string(name) + " is "
					+ std::to_string(a) + " in the SOLEDGE3X reference file, "
					"but " + std::to_string(b) + " in the mesh file. They "
					"should be from the same case.");
			}
		}
	}

	/**
	* Reads the number of toroidal planes (Nphi) of a SOLEDGE3X plasma file.
	*
	* Every zone must have the same Nphi, otherwise the run is stopped.
	*
	* @param plasma_path Path to a plasma file.
	* @param nzones Number of zones.
	* @return Nphi.
	*/
	int read_nphi(const std::string& plasma_path, const int nzones)
	{
		const Hdf5::File file(plasma_path);
		int nphi {};
		for (int zone {}; zone < nzones; ++zone)
		{
			const int n {Hdf5::read_scalar<int>(file, zone_group(zone)
				+ "/Nphi")};
			if (zone == 0) nphi = n;
			else if (n != nphi)
			{
				Errors::fatal_error("In " + plasma_path + ", zone "
					+ std::to_string(zone + 1) + " has Nphi = "
					+ std::to_string(n) + ", but zone 1 has "
					+ std::to_string(nphi) + ".");
			}
		}
		if (nphi < 1)
		{
			Errors::fatal_error("In " + plasma_path + ", Nphi is "
				+ std::to_string(nphi) + ".");
		}
		return nphi;
	}

	/**
	* Reads a SOLEDGE3X mesh file and builds its zones.
	*
	* The number of zones comes from NZones and must match the neighbors
	* table. Each zone gets x = psi cell index in [0, Npsi], y = theta cell
	* index in [0, Ntheta] and z = phi with nphi cells between the zone's
	* phicorners, and ghost cells of width ghost. Its faces come from the
	* neighbors table, whose rows are psi-, psi+, theta-, theta+, phi- and
	* phi+: a zone number (from 1) connects to that zone, -2 is the core,
	* and -1 (target) and -3 (outer boundary) absorb particles. The zones
	* are checked with Zones::Zones::validate.
	*
	* Only axisymmetric meshes (Nphi = 1 in the mesh file) are handled. Any
	* missing dataset, wrong shape or unknown neighbor code stops the run.
	*
	* @param mesh_path Path to the mesh file.
	* @param nphi Number of toroidal planes of the plasma (see read_nphi).
	* @return The mesh, still normalized, with its zones.
	*/
	Mesh read_mesh(const std::string& mesh_path, const int nphi)
	{
		const Hdf5::File file(mesh_path);
		Mesh mesh {};
		mesh.R0 = Hdf5::read_scalar<double>(file, "/R0");
		mesh.a0 = Hdf5::read_scalar<double>(file, "/a0");
		mesh.B0 = Hdf5::read_scalar<double>(file, "/B0");

		const int nzones {Hdf5::read_scalar<int>(file, "/NZones")};
		if (nzones < 1)
		{
			Errors::fatal_error("In " + mesh_path + ", NZones is "
				+ std::to_string(nzones) + ".");
		}
		const auto neighbors {Hdf5::read_dataset<int>(file, "/neighbors")};
		check_shape(neighbors.shape, {6, static_cast<hsize_t>(nzones)},
			"/neighbors", mesh_path);

		mesh.zones = Zones::Zones(ghost);
		for (int zone {}; zone < nzones; ++zone)
		{
			const std::string g {zone_group(zone)};
			const int npsi {Hdf5::read_scalar<int>(file, g + "/Npsi")};
			const int ntheta {Hdf5::read_scalar<int>(file, g + "/Ntheta")};
			const int nphi_mesh {Hdf5::read_scalar<int>(file, g + "/Nphi")};
			if (npsi < 1 || ntheta < 1)
			{
				Errors::fatal_error("In " + mesh_path + ", zone "
					+ std::to_string(zone + 1) + " has Npsi = "
					+ std::to_string(npsi) + " and Ntheta = "
					+ std::to_string(ntheta) + ".");
			}
			if (nphi_mesh != 1)
			{
				Errors::fatal_error("In " + mesh_path + ", zone "
					+ std::to_string(zone + 1) + " has Nphi = "
					+ std::to_string(nphi_mesh) + ". Only axisymmetric "
					"meshes (Nphi = 1) are handled.");
			}

			const auto phic {read_checked(file, g + "/phicorners", {2})};
			if (!(phic[1] > phic[0]))
			{
				Errors::fatal_error("In " + mesh_path + ", zone "
					+ std::to_string(zone + 1) + "'s phicorners are not "
					"increasing.");
			}

			const hsize_t np {static_cast<hsize_t>(npsi)};
			const hsize_t nt {static_cast<hsize_t>(ntheta)};
			mesh.R_corners.push_back(read_checked(file, g + "/Rcorners",
				{np + 1, nt + 1, 1}));
			mesh.Z_corners.push_back(read_checked(file, g + "/Zcorners",
				{np + 1, nt + 1, 1}));
			mesh.Br.push_back(read_checked(file, g + "/Br", {np, nt, 1}));
			mesh.Bz.push_back(read_checked(file, g + "/Bz", {np, nt, 1}));
			mesh.Bphi.push_back(read_checked(file, g + "/Bphi",
				{np, nt, 1}));
			mesh.chi.push_back(read_checked(file, g + "/chi", {np, nt, 1}));

			auto [x, gx] = uniform_grid(npsi, 0.0, npsi);
			auto [y, gy] = uniform_grid(ntheta, 0.0, ntheta);
			auto [z, gz] = uniform_grid(nphi, phic[0], phic[1]);
			mesh.zones.add_zone(x, y, z, gx, gy, gz);
		}

		// Faces from the neighbors table, stored as (6, NZones)
		for (int zone {}; zone < nzones; ++zone)
		{
			for (int face {}; face < Zones::num_faces; ++face)
			{
				const int code {neighbors.data[face * nzones + zone]};
				int face_code {};
				if (code > 0 && code <= nzones) face_code = code - 1;
				else if (code == -2) face_code = Zones::face_core;
				else if (code == -1 || code == -3)
					face_code = Zones::face_absorbing;
				else
				{
					Errors::fatal_error("In " + mesh_path + ", zone "
						+ std::to_string(zone + 1) + " has neighbor code "
						+ std::to_string(code) + ", which isn't a zone "
						"(1 to " + std::to_string(nzones) + "), -1, -2 or "
						"-3.");
				}
				mesh.zones.set_face(zone, static_cast<Zones::Face>(face),
					face_code);
			}
		}
		mesh.zones.validate();

		return mesh;
	}

	/**
	* Derivatives of a field with respect to the computational coordinates
	* x, y and z at an interior cell.
	*
	* A central difference is used where the cells on both sides exist,
	* including ghost cells across faces shared with another zone (or
	* periodic), so ghost cells must be filled. At a face with nothing
	* across it (absorbing or core), the difference is one-sided. Cell
	* spacings come from Zones::Zones::center, so they're right across
	* zones too.
	*
	* @param zones The zones.
	* @param f Field, with filled ghost cells.
	* @param zone Zone of the cell.
	* @param t Index along f's leading dimension.
	* @param i Cell index in x.
	* @param j Cell index in y.
	* @param k Cell index in z.
	* @return {df/dx, df/dy, df/dz}.
	*/
	std::array<double, 3> derivatives(const Zones::Zones& zones,
		const Vectors::ZonedVector<BkgFPType>& f, const int zone,
		const int t, const int i, const int j, const int k)
	{
		const std::array<int, 3> idx {i, j, k};
		const bool ghosts {zones.get_ghost() > 0};
		std::array<double, 3> out {};
		for (int d {}; d < 3; ++d)
		{
			const int n {zones.get_n(zone, d)};
			const bool has_lo {idx[d] > 0 || (ghosts
				&& zones.get_face(zone, static_cast<Zones::Face>(2 * d)) >= 0)};
			const bool has_hi {idx[d] < n - 1 || (ghosts
				&& zones.get_face(zone, static_cast<Zones::Face>(2 * d + 1))
				>= 0)};
			const int lo {has_lo ? idx[d] - 1 : idx[d]};
			const int hi {has_hi ? idx[d] + 1 : idx[d]};
			if (lo == hi) continue;

			std::array<int, 3> a {idx};
			std::array<int, 3> b {idx};
			a[d] = lo;
			b[d] = hi;
			out[d] = (f(zone, t, b[0], b[1], b[2]) - f(zone, t, a[0], a[1],
				a[2])) / (zones.center(zone, d, hi) - zones.center(zone, d, lo));
		}
		return out;
	}

	/**
	* Cartesian gradient of a field at an interior cell, from its
	* derivatives (see derivatives) and the background's reciprocal basis:
	* grad f = df/dx grad x + df/dy grad y + df/dz grad z.
	*
	* @param bkg Background holding the zones and the reciprocal basis.
	* @param f Field, with filled ghost cells.
	* @param zone Zone of the cell.
	* @param t Index along f's leading dimension.
	* @param i Cell index in x.
	* @param j Cell index in y.
	* @param k Cell index in z.
	* @return {df/dX, df/dY, df/dZ}.
	*/
	std::array<double, 3> gradient(const Background::Background& bkg,
		const Vectors::ZonedVector<BkgFPType>& f, const int zone,
		const int t, const int i, const int j, const int k)
	{
		const auto df {derivatives(bkg.get_zones(), f, zone, t, i, j, k)};
		const std::array<std::array<const Vectors::ZonedVector<BkgFPType>*,
			3>, 3> recp {{
			{&bkg.get_dxdX(), &bkg.get_dxdY(), &bkg.get_dxdZ()},
			{&bkg.get_dydX(), &bkg.get_dydY(), &bkg.get_dydZ()},
			{&bkg.get_dzdX(), &bkg.get_dzdY(), &bkg.get_dzdZ()}}};

		std::array<double, 3> grad {};
		for (int d {}; d < 3; ++d)
		for (int c {}; c < 3; ++c)
		{
			grad[c] += df[d] * (*recp[d][c])(zone, i, j, k);
		}
		return grad;
	}

	/**
	* Calculates the geometry of every cell from the mesh corners and moves
	* it into the background.
	*
	* The background's zones must already be set to mesh.zones. In each
	* toroidal plane, a cell's center (R, Z) is the bilinear center of its
	* four corners, in meters. The tangent basis vectors along x and y are
	* the differences between the midpoints of opposite cell faces, and
	* along z it's R times the unit vector in phi. From these come the
	* Jacobian J = e_x . (e_y x e_z) and the reciprocal basis, e.g.,
	* grad x = (e_y x e_z) / J.
	*
	* Moved into the background are the cell centers X, Y and Z, the
	* reciprocal basis (dxdX, ..., dzdZ) and |J|, since psi, theta and phi
	* may be a left-handed set but J is used as a cell volume, all with
	* their ghost cells filled, and the cell corners grid_X, grid_Y and
	* grid_Z. The tangent basis and metric aren't kept, to save memory.
	*
	* @param mesh The mesh.
	* @param bkg Background to move the geometry into.
	*/
	void calc_geometry(const Mesh& mesh, Background::Background& bkg)
	{
		const Zones::Zones& zones {bkg.get_zones()};
		auto X {zoned(zones, 1)};
		auto Y {zoned(zones, 1)};
		auto Z {zoned(zones, 1)};
		auto J {zoned(zones, 1)};
		std::array<Vectors::ZonedVector<BkgFPType>, 9> recp {};
		for (auto& r : recp) r = zoned(zones, 1);

		// Cell corners have one more value in each dimension and no ghost
		// cells
		auto corners = [&]()
		{
			std::vector<int> nx {}, ny {}, nz {};
			for (int zone {}; zone < zones.get_nzones(); ++zone)
			{
				nx.push_back(zones.get_nx(zone) + 1);
				ny.push_back(zones.get_ny(zone) + 1);
				nz.push_back(zones.get_nz(zone) + 1);
			}
			return Vectors::ZonedVector<BkgFPType>(1, nx, ny, nz, 0);
		};
		auto grid_X {corners()};
		auto grid_Y {corners()};
		auto grid_Z {corners()};

		for (int zone {}; zone < zones.get_nzones(); ++zone)
		{
			const int npsi {zones.get_nx(zone)};
			const int ntheta {zones.get_ny(zone)};
			const int nphi {zones.get_nz(zone)};
			const auto& Rc {mesh.R_corners[zone]};
			const auto& Zc {mesh.Z_corners[zone]};

			// Corner (i, j) in meters
			auto R = [&](int i, int j)
				{return mesh.R0 + mesh.a0 * Rc[i * (ntheta + 1) + j];};
			auto Zm = [&](int i, int j)
				{return mesh.a0 * Zc[i * (ntheta + 1) + j];};

			for (int i {}; i < npsi; ++i)
			for (int j {}; j < ntheta; ++j)
			{
				const double R_cen {0.25 * (R(i, j) + R(i + 1, j)
					+ R(i, j + 1) + R(i + 1, j + 1))};
				const double Z_cen {0.25 * (Zm(i, j) + Zm(i + 1, j)
					+ Zm(i, j + 1) + Zm(i + 1, j + 1))};
				const double dRdx {0.5 * (R(i + 1, j) + R(i + 1, j + 1)
					- R(i, j) - R(i, j + 1))};
				const double dZdx {0.5 * (Zm(i + 1, j) + Zm(i + 1, j + 1)
					- Zm(i, j) - Zm(i, j + 1))};
				const double dRdy {0.5 * (R(i, j + 1) + R(i + 1, j + 1)
					- R(i, j) - R(i + 1, j))};
				const double dZdy {0.5 * (Zm(i, j + 1) + Zm(i + 1, j + 1)
					- Zm(i, j) - Zm(i + 1, j))};

				for (int k {}; k < nphi; ++k)
				{
					const double phi {zones.get_z(zone)[k]};
					const Vec3 e_x {combine(dRdx, R_hat(phi), dZdx, Z_hat)};
					const Vec3 e_y {combine(dRdy, R_hat(phi), dZdy, Z_hat)};
					const Vec3 e_z {combine(R_cen, phi_hat(phi), 0.0,
						Z_hat)};
					const double jac {dot(e_x, cross(e_y, e_z))};
					const std::array<Vec3, 3> r {cross(e_y, e_z),
						cross(e_z, e_x), cross(e_x, e_y)};

					X(zone, i, j, k) = R_cen * std::cos(phi);
					Y(zone, i, j, k) = phi_sign * R_cen * std::sin(phi);
					Z(zone, i, j, k) = Z_cen;
					J(zone, i, j, k) = std::abs(jac);
					for (int d {}; d < 3; ++d)
					for (int c {}; c < 3; ++c)
					{
						recp[3 * d + c](zone, i, j, k) = r[d][c] / jac;
					}
				}
			}

			for (int i {}; i <= npsi; ++i)
			for (int j {}; j <= ntheta; ++j)
			for (int k {}; k <= nphi; ++k)
			{
				const double phi {zones.get_grid_z(zone)[k]};
				grid_X(zone, i, j, k) = R(i, j) * std::cos(phi);
				grid_Y(zone, i, j, k) = phi_sign * R(i, j) * std::sin(phi);
				grid_Z(zone, i, j, k) = Zm(i, j);
			}
		}

		for (auto* v : {&X, &Y, &Z, &J}) zones.fill_ghosts(*v);
		for (auto& r : recp) zones.fill_ghosts(r);

		bkg.move_into_X(X);
		bkg.move_into_Y(Y);
		bkg.move_into_Z(Z);
		bkg.move_into_J(J);
		bkg.move_into_grid_X(grid_X);
		bkg.move_into_grid_Y(grid_Y);
		bkg.move_into_grid_Z(grid_Z);
		bkg.move_into_dxdX(recp[0]);
		bkg.move_into_dxdY(recp[1]);
		bkg.move_into_dxdZ(recp[2]);
		bkg.move_into_dydX(recp[3]);
		bkg.move_into_dydY(recp[4]);
		bkg.move_into_dydZ(recp[5]);
		bkg.move_into_dzdX(recp[6]);
		bkg.move_into_dzdY(recp[7]);
		bkg.move_into_dzdZ(recp[8]);
	}

	/**
	* Calculates the magnetic field and the gradient of its magnitude, and
	* moves them into the background.
	*
	* The mesh's (Br, Bphi, Bz) times B0 are rotated into Cartesian
	* components in each toroidal plane. The field is axisymmetric and
	* doesn't change in time, so every plane and every one of the nframes
	* time frames get the same values. grad |B| is calculated with gradient,
	* so calc_geometry must be called first. Moved into the background are
	* bX, bY, bZ, bmag, gradbX, gradbY and gradbZ, with ghost cells filled.
	*
	* @param mesh The mesh.
	* @param nframes Number of time frames.
	* @param bkg Background to move the field into, with its geometry set.
	*/
	void calc_bfield(const Mesh& mesh, const int nframes,
		Background::Background& bkg)
	{
		if (bkg.get_dxdX().get_data().empty())
		{
			Errors::fatal_error("S3X::calc_bfield was called before "
				"S3X::calc_geometry.");
		}

		const Zones::Zones& zones {bkg.get_zones()};
		auto bX {zoned(zones, nframes)};
		auto bY {zoned(zones, nframes)};
		auto bZ {zoned(zones, nframes)};
		auto bmag {zoned(zones, nframes)};
		auto gradbX {zoned(zones, nframes)};
		auto gradbY {zoned(zones, nframes)};
		auto gradbZ {zoned(zones, nframes)};

		for (int zone {}; zone < zones.get_nzones(); ++zone)
		{
			const int ntheta {zones.get_ny(zone)};
			for (int i {}; i < zones.get_nx(zone); ++i)
			for (int j {}; j < ntheta; ++j)
			{
				const int c {i * ntheta + j};
				const double BR {mesh.B0 * mesh.Br[zone][c]};
				const double BZ {mesh.B0 * mesh.Bz[zone][c]};
				const double Bphi {mesh.B0 * mesh.Bphi[zone][c]};
				for (int k {}; k < zones.get_nz(zone); ++k)
				{
					const double phi {zones.get_z(zone)[k]};
					const Vec3 B {combine(BR, R_hat(phi), Bphi,
						phi_hat(phi))};
					for (int t {}; t < nframes; ++t)
					{
						bX(zone, t, i, j, k) = B[0];
						bY(zone, t, i, j, k) = B[1];
						bZ(zone, t, i, j, k) = BZ;
						bmag(zone, t, i, j, k) = std::sqrt(B[0] * B[0]
							+ B[1] * B[1] + BZ * BZ);
					}
				}
			}
		}
		for (auto* v : {&bX, &bY, &bZ, &bmag}) zones.fill_ghosts(*v);

		// Gradient of |B|, the same in every time frame
		for (int zone {}; zone < zones.get_nzones(); ++zone)
		for (int i {}; i < zones.get_nx(zone); ++i)
		for (int j {}; j < zones.get_ny(zone); ++j)
		for (int k {}; k < zones.get_nz(zone); ++k)
		{
			const auto grad {gradient(bkg, bmag, zone, 0, i, j, k)};
			for (int t {}; t < nframes; ++t)
			{
				gradbX(zone, t, i, j, k) = grad[0];
				gradbY(zone, t, i, j, k) = grad[1];
				gradbZ(zone, t, i, j, k) = grad[2];
			}
		}
		for (auto* v : {&gradbX, &gradbY, &gradbZ}) zones.fill_ghosts(*v);

		bkg.move_into_bX(bX);
		bkg.move_into_bY(bY);
		bkg.move_into_bZ(bZ);
		bkg.move_into_bmag(bmag);
		bkg.move_into_gradbX(gradbX);
		bkg.move_into_gradbY(gradbY);
		bkg.move_into_gradbZ(gradbZ);
	}

	/**
	* Moves the wall mask chi into the background, the same in every
	* toroidal plane, with ghost cells filled. chi is 1 in cells inside the
	* wall and 0 in the plasma.
	*
	* @param mesh The mesh.
	* @param bkg Background to move chi into, with its zones set.
	*/
	void load_chi(const Mesh& mesh, Background::Background& bkg)
	{
		const Zones::Zones& zones {bkg.get_zones()};
		auto chi {zoned(zones, 1)};
		for (int zone {}; zone < zones.get_nzones(); ++zone)
		{
			const int ntheta {zones.get_ny(zone)};
			for (int i {}; i < zones.get_nx(zone); ++i)
			for (int j {}; j < ntheta; ++j)
			for (int k {}; k < zones.get_nz(zone); ++k)
			{
				chi(zone, i, j, k) = mesh.chi[zone][i * ntheta + j];
			}
		}
		zones.fill_ghosts(chi);
		bkg.move_into_chi(chi);
	}
}
