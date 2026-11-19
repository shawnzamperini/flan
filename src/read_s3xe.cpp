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
* for a single toroidal plane (the mesh is axisymmetric). The mesh is
* converted to SI units (to_axisym), and the geometry, magnetic field,
* electric field and flows are computed from it by the code shared with
* other axisymmetric multi-zone backgrounds (axisym_mesh.cpp).
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "axisym_mesh.h"
#include "background.h"
#include "constants.h"
#include "errors.h"
#include "flan_types.h"
#include "hdf5_utils.h"
#include "options.h"
#include "read_s3xe.h"
#include "vectors.h"
#include "zones.h"

namespace S3X
{
	namespace
	{
		// Name of a zone's group in the SOLEDGE3X files, numbered from 1
		std::string zone_group(const int zone)
		{
			return "/zone" + std::to_string(zone + 1);
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

		/// A phi face used to fit C: SOLEDGE3X's particle flux through it
		/// (fluxn) and the flux density along B times the face's geometry
		/// (G b_phi A_pol / a0^2), which C turns into fluxn.
		struct FitPoint
		{
			double fluxn {};
			double G_geom {};
		};

		/**
		* Fits the constant C in fluxn = C G b_phi A_pol / a0^2, which turns a
		* parallel particle flux density G into SOLEDGE3X's face flux (see
		* calc_q_par_i). C is the median of fluxn / (G b_phi A_pol / a0^2)
		* over the faces whose |G b_phi A_pol / a0^2| is at least a tenth of
		* the largest, so faces with almost no parallel flux, where turbulent
		* fluxes and round-off dominate the ratio, are left out. The median
		* also ignores the few faces where they still do.
		*
		* @param points Faces in the plasma.
		* @param path Plasma file, for the error message.
		* @return C.
		*/
		double fit_face_flux_constant(const std::vector<FitPoint>& points,
			const std::string& path)
		{
			double largest {};
			for (const auto& p : points)
				largest = std::max(largest, std::abs(p.G_geom));

			std::vector<double> ratios {};
			for (const auto& p : points)
			{
				if (largest > 0.0 && std::abs(p.G_geom) >= 0.1 * largest)
					ratios.push_back(p.fluxn / p.G_geom);
			}
			if (ratios.size() < 10)
			{
				Errors::fatal_error("Can't fit the constant relating "
					"SOLEDGE3X's phi face fluxes to G in " + path + ": only "
					+ std::to_string(ratios.size()) + " faces have a "
					"parallel flux.");
			}
			const auto mid {ratios.begin() + ratios.size() / 2};
			std::nth_element(ratios.begin(), mid, ratios.end());
			return *mid;
		}
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
	* The number of zones comes from the file's NZones. Every zone must have
	* the same Nphi, otherwise the run is stopped.
	*
	* @param plasma_path Path to a plasma file.
	* @return Nphi.
	*/
	int read_nphi(const std::string& plasma_path)
	{
		const Hdf5::File file(plasma_path);
		const int nzones {Hdf5::read_scalar<int>(file, "/NZones")};
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
	* -1 (target) is a wall face and -3 (outer boundary) is an absorbing
	* face. The zones are checked with Zones::Zones::validate.
	*
	* The wall outlines in the /walls group are read too, if the file has
	* one: Nwalls outlines wall1, wall2, ..., each with R and Z in meters
	* (not normalized like the corners). Their type isn't used. A file
	* without /walls has no outlines.
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

			auto [x, gx] = Axisym::uniform_grid(npsi, 0.0, npsi);
			auto [y, gy] = Axisym::uniform_grid(ntheta, 0.0, ntheta);
			auto [z, gz] = Axisym::uniform_grid(nphi, phic[0], phic[1]);
			mesh.zones.add_zone(x, y, z, gx, gy, gz);

			// z is the toroidal angle, so Cartesian vectors turn across
			// the phi faces of a zone narrower than the whole torus
			mesh.zones.set_toroidal(zone);
		}

		// Faces from the neighbors table, stored as (6, NZones). SOLEDGE3X
		// codes: > 0 a zone (from 1), -1 a divertor target (a wall face in
		// Flan), -2 the core and -3 the outer psi boundary (an absorbing
		// face in Flan).
		for (int zone {}; zone < nzones; ++zone)
		{
			for (int face {}; face < Zones::num_faces; ++face)
			{
				const int code {neighbors.data[face * nzones + zone]};
				int face_code {};
				if (code > 0 && code <= nzones) face_code = code - 1;
				else if (code == -2) face_code = Zones::face_core;
				else if (code == -1) face_code = Zones::face_wall;
				else if (code == -3) face_code = Zones::face_absorbing;
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

		// Wall outlines, if there are any. They're already in meters.
		if (file.exists("/walls/Nwalls"))
		{
			const int nwalls {static_cast<int>(Hdf5::read_scalar<double>(
				file, "/walls/Nwalls"))};
			for (int w {1}; w <= nwalls; ++w)
			{
				const std::string g {"/walls/wall" + std::to_string(w)};
				Background::WallPolygon wall {};
				wall.R = Hdf5::read_dataset<double>(file, g + "/R").data;
				wall.Z = Hdf5::read_dataset<double>(file, g + "/Z").data;
				if (wall.R.size() != wall.Z.size() || wall.R.size() < 2)
				{
					Errors::fatal_error("In " + mesh_path + ", " + g
						+ " has " + std::to_string(wall.R.size()) + " R and "
						+ std::to_string(wall.Z.size()) + " Z values. They "
						"need to be the same, and at least 2.");
				}
				mesh.walls.push_back(std::move(wall));
			}
		}

		return mesh;
	}

	/**
	* Converts the mesh's corners and magnetic field to SI units for the
	* shared geometry code: R = R0 + a0 Rcorners, Z = a0 Zcorners (m) and
	* B = B0 (Br, Bz, Bphi) (T), with SOLEDGE3X's handedness (phi_sign).
	* The array layouts are already the ones Axisym::Mesh uses.
	*
	* @param mesh The mesh, as read by read_mesh.
	* @return The mesh in SI units.
	*/
	Axisym::Mesh to_axisym(const Mesh& mesh)
	{
		Axisym::Mesh out {};
		out.phi_sign = phi_sign;
		// factor * value. R has its own loop for the offset, since adding a
		// zero offset to Z and B would turn -0 into +0.
		auto scale = [](const std::vector<std::vector<double>>& in,
			const double factor)
		{
			std::vector<std::vector<double>> v {in};
			for (auto& zone : v)
				for (auto& val : zone) val = factor * val;
			return v;
		};
		out.R_corners = mesh.R_corners;
		for (auto& zone : out.R_corners)
			for (auto& val : zone) val = mesh.R0 + mesh.a0 * val;
		out.Z_corners = scale(mesh.Z_corners, mesh.a0);
		out.Br = scale(mesh.Br, mesh.B0);
		out.Bz = scale(mesh.Bz, mesh.B0);
		out.Bphi = scale(mesh.Bphi, mesh.B0);
		return out;
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
		auto chi {Axisym::zoned(zones, 1)};
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

	/**
	* Whether a SOLEDGE3X plasma file has a species' fluxes on the phi
	* faces (e.g., spec1/fluxn/phi and spec1/fluxE/phi), which only runs set
	* up to save face fluxes have. Only zone 1 is checked; if it has them,
	* every zone must, or reading them stops the run.
	*
	* @param file Plasma file.
	* @param species The species' group: spec0 (electrons) or spec1 (ions).
	* @return True if the file has them.
	*/
	bool has_phi_face_fluxes(const Hdf5::File& file,
		const std::string& species)
	{
		return file.exists(zone_group(0) + "/" + species + "/fluxn/phi")
			&& file.exists(zone_group(0) + "/" + species + "/fluxE/phi");
	}

	/**
	* Calculates a species' parallel conductive heat flux of one frame from
	* SOLEDGE3X's fluxes of that species through the phi faces. It is the
	* same for electrons (spec0) and ions (spec1); only the group and mass
	* differ.
	*
	* What is read, and why phi faces: SOLEDGE3X saves, for every face of
	* every cell, the particle flux fluxn and the total energy flux fluxE
	* crossing it. Through a phi face (a face in the poloidal plane), B is
	* almost normal (b_phi ~ 1), so these are almost entirely the parallel
	* fluxes. Through theta faces, B is at a shallow angle (b . n ~ 0.1-0.3)
	* and the turbulent and drift fluxes crossing them rival the parallel
	* part, so they can't be used.
	*
	* Total vs. conductive: fluxE is the total energy flux, carried by the
	* flow (the enthalpy 5/2 T and the kinetic energy 1/2 m u^2 per
	* particle) and by conduction. Only the conductive heat flux q_par is
	* kept: the thermal force model samples the species' distribution in
	* the frame moving with the flow (the flow itself is already the drift of
	* the Maxwellian), and in that frame only the conductive part remains
	* as a skew of the distribution. Per face,
	*
	*   q_par = (fluxE - (5/2 T + 1/2 m u^2) fluxn) / (C b_phi A_pol / a0^2)
	*
	* where, in SOLEDGE3X's normalized units:
	*   - T, u: the species' temperature (T0) and parallel velocity G / n
	*     (c0) on the face, from the average of n, T and G in the two cells
	*     on either side. m is its mass in amu, so 1/2 m u^2 is in T0 as
	*     well, since c0 = sqrt(e T0 / 1 amu). (For electrons it's
	*     negligible next to 5/2 T.)
	*   - b_phi = Bphi / |B|, the component of the unit vector along B
	*     normal to the face.
	*   - A_pol / a0^2: the face's area, which is the cell's area in the
	*     poloidal plane, from the shoelace formula on the cell's four
	*     normalized corners (the R0 offset doesn't change an area).
	*   - C: the constant that turns a flux density along B into what
	*     SOLEDGE3X saves as a face flux, fluxn = C G b_phi A_pol / a0^2.
	*     Its definition isn't documented, so it's fitted here (see
	*     fit_face_flux_constant), for each species separately as a check:
	*     it's geometric, so they should agree. For the DIII-D runs it is
	*     close to -2 pi c0 / a0, which is printed for comparison.
	*
	* Dividing by C b_phi A_pol / a0^2 turns the face flux back into a
	* flux density along B. Since C is fitted against G, which is positive
	* along B, q_par is positive along B too, like u_par.
	*
	* The cell-centered value is the average of the cell's two phi faces,
	* in W/m^2 (q_par n0 c0 T0 e). Inside the wall (chi = 1), SOLEDGE3X's
	* plasma is penalized and q_par there means nothing physically.
	*
	* @param file Plasma file of this frame.
	* @param ref Reference values.
	* @param mesh The mesh, for the corners and B.
	* @param zones The zones.
	* @param species The species' group: spec0 (electrons) or spec1 (ions).
	* @param mass_amu The species' mass (amu).
	* @param t Frame index.
	* @param q_par Where q_par goes, for frame t (ghost cells not filled).
	* @param strict If true, a fitted C more than 5% from -2 pi c0 / a0
	*        stops the run (the thermal force uses q_par); otherwise it's a
	*        warning.
	* @return The fitted C, for printing.
	*/
	double calc_q_par(const Hdf5::File& file, const Reference& ref,
		const Mesh& mesh, const Zones::Zones& zones,
		const std::string& species, const double mass_amu, const int t,
		Vectors::ZonedVector<BkgFPType>& q_par, const bool strict)
	{
		const std::string sp {"/" + species};
		const int g {zones.get_ghost()};
		const int nzones {zones.get_nzones()};

		// Per zone and face (psi i, theta j, phi face k = 0..Nphi, face k
		// between cells k - 1 and k): the conductive part of the energy
		// flux in SOLEDGE3X's units (fluxE - h fluxn), b_phi A_pol / a0^2,
		// and what the fit of C needs.
		std::vector<std::vector<double>> conductive(nzones);
		std::vector<std::vector<double>> geom(nzones);
		std::vector<FitPoint> fit_points {};

		for (int zone {}; zone < nzones; ++zone)
		{
			const std::string zg {zone_group(zone)};
			const int nx {zones.get_nx(zone)};
			const int ny {zones.get_ny(zone)};
			const int nz {zones.get_nz(zone)};
			const std::vector<hsize_t> cells {
				static_cast<hsize_t>(nx + 2 * g),
				static_cast<hsize_t>(ny + 2 * g),
				static_cast<hsize_t>(nz + 2 * g)};
			const std::vector<hsize_t> faces {static_cast<hsize_t>(nx),
				static_cast<hsize_t>(ny), static_cast<hsize_t>(nz + 1)};
			const auto n {read_checked(file, zg + sp + "/n", cells)};
			const auto T {read_checked(file, zg + sp + "/T", cells)};
			const auto G {read_checked(file, zg + sp + "/G", cells)};
			const auto fluxn {read_checked(file, zg + sp + "/fluxn/phi",
				faces)};
			const auto fluxE {read_checked(file, zg + sp + "/fluxE/phi",
				faces)};

			// Cell (i, j, k) in the plasma file's arrays, which have ghost
			// cells, so k = -1 is the ghost cell before the first plane
			auto cell = [&](const int i, const int j, const int k)
			{
				return (static_cast<std::size_t>(i + g) * (ny + 2 * g)
					+ (j + g)) * (nz + 2 * g) + (k + g);
			};

			conductive[zone].resize(fluxn.size());
			geom[zone].resize(fluxn.size());
			const auto& Rc {mesh.R_corners[zone]};
			const auto& Zc {mesh.Z_corners[zone]};
			for (int i {}; i < nx; ++i)
			for (int j {}; j < ny; ++j)
			{
				// Poloidal area of the cell (normalized by a0^2), with the
				// shoelace formula on its corners in order around it
				const std::array<std::size_t, 4> c {
					static_cast<std::size_t>(i * (ny + 1) + j),
					static_cast<std::size_t>((i + 1) * (ny + 1) + j),
					static_cast<std::size_t>((i + 1) * (ny + 1) + j + 1),
					static_cast<std::size_t>(i * (ny + 1) + j + 1)};
				double twice_area {};
				for (int p {}; p < 4; ++p)
				{
					const std::size_t a {c[p]};
					const std::size_t b {c[(p + 1) % 4]};
					twice_area += Rc[a] * Zc[b] - Rc[b] * Zc[a];
				}
				const double A_pol {0.5 * std::abs(twice_area)};

				const int m {i * ny + j};
				const double Bmag {std::sqrt(mesh.Br[zone][m]
					* mesh.Br[zone][m] + mesh.Bz[zone][m] * mesh.Bz[zone][m]
					+ mesh.Bphi[zone][m] * mesh.Bphi[zone][m])};
				const double b_phi {Bmag > 0.0 ? mesh.Bphi[zone][m] / Bmag
					: 0.0};
				const bool plasma {mesh.chi[zone][m] < 0.5};

				for (int k {}; k <= nz; ++k)
				{
					const std::size_t lo {cell(i, j, k - 1)};
					const std::size_t hi {cell(i, j, k)};
					const double n_face {0.5 * (n[lo] + n[hi])};
					const double T_face {0.5 * (T[lo] + T[hi])};
					const double G_face {0.5 * (G[lo] + G[hi])};
					const double u_face {n_face > 0.0 ? G_face / n_face
						: 0.0};

					// Energy carried by the flow per particle: enthalpy
					// 5/2 T plus kinetic energy 1/2 m u^2 (both in T0)
					const double h {2.5 * T_face
						+ 0.5 * mass_amu * u_face * u_face};

					const std::size_t f {(static_cast<std::size_t>(i) * ny
						+ j) * (nz + 1) + k};
					conductive[zone][f] = fluxE[f] - h * fluxn[f];
					geom[zone][f] = b_phi * A_pol;
					if (plasma)
					{
						fit_points.push_back({fluxn[f], G_face
							* geom[zone][f]});
					}
				}
			}
		}

		// Unnormalize to W/m^2
		const double C {fit_face_flux_constant(fit_points, file.path())};
		const double q0 {ref.n0 * ref.c0 * ref.T0 * Constants::ev_to_j};
		if (t == 0)
		{
			std::cout << "  Face flux constant C = " << C << " from "
				<< species << " (-2 pi c0 / a0 = "
				<< -2.0 * Constants::pi * ref.c0 / mesh.a0 << ")\n";
		}
		const double C_expect {2.0 * Constants::pi * ref.c0 / mesh.a0};
		if (std::abs(std::abs(C) / C_expect - 1.0) > 0.05 && strict)
		{
			Errors::fatal_error("In " + file.path() + ", the fitted face "
				"flux constant C = " + std::to_string(C) + " from " + species
				+ " is more than 5% from -2 pi c0 / a0 = "
				+ std::to_string(-C_expect) + ", so SOLEDGE3X's face fluxes "
				"may not be defined as expected. The thermal force "
				"(imp_thermal_force = \"on\") needs them; set it to \"off\" "
				"to run without it.");
		}
		if (std::abs(std::abs(C) / C_expect - 1.0) > 0.05)
		{
			std::cout << "Warning! In " << file.path() << ", the fitted "
				"face flux constant C = " << C << " from " << species
				<< " is more than 5% from -2 pi c0 / a0 = " << -C_expect
				<< ". The heat flux still uses the fitted C, but check that "
				"SOLEDGE3X's face fluxes are defined as expected.\n";
		}

		for (int zone {}; zone < nzones; ++zone)
		{
			const int ny {zones.get_ny(zone)};
			const int nz {zones.get_nz(zone)};
			for (int i {}; i < zones.get_nx(zone); ++i)
			for (int j {}; j < ny; ++j)
			for (int k {}; k < nz; ++k)
			{
				// q_par on the cell's two phi faces, then their average
				double sum {};
				for (const int kf : {k, k + 1})
				{
					const std::size_t f {(static_cast<std::size_t>(i) * ny
						+ j) * (nz + 1) + kf};
					const double denom {C * geom[zone][f]};
					sum += denom != 0.0 ? conductive[zone][f] / denom : 0.0;
				}
				q_par(zone, t, i, j, k) = 0.5 * sum * q0;
			}
		}
		return C;
	}

	/**
	* Reads the plasma of every frame and moves it into the background.
	*
	* The background must already have its zones, geometry and magnetic
	* field. For each frame, the electron density and temperature (spec0),
	* ion temperature (spec1) and potential (PHI) are read and unnormalized
	* with the reference values: ne = n n0, Te = T T0, Ti = T T0 and
	* vp = PHI phi0. The parallel ion flow is u_par = G / n c0, from the
	* ion's G and n. The time is tsimu tau0, and the times must increase.
	* The masses of the electrons and ions come from plasma_composition in
	* the first frame.
	*
	* The plasma files already have ghost cells. How they differ from
	* fill_ghosts is printed for the first frame (see report_ghosts), and
	* then they are replaced with fill_ghosts, so the plasma uses the same
	* ghost cells as the geometry. Then E and the flows are calculated with
	* calc_efield_and_flows.
	*
	* For each species, if the first frame has its fluxes through the phi
	* faces (runs set up to save face fluxes), its parallel conductive heat
	* flux (q_par_e for electrons, q_par_i for ions) is calculated from them
	* (see calc_q_par) and moved into the background too, and every frame
	* must have them. Otherwise it is left empty and a message says so.
	*
	* Any missing dataset, wrong shape, or a number of zones or ghost cells
	* that doesn't match the mesh stops the run.
	*
	* @param frame_paths Paths of the plasma files, in time order.
	* @param ref Reference values.
	* @param mesh The mesh, for the heat flux.
	* @param bkg Background to move the plasma into.
	* @param strict If true, a face flux constant far from what's expected
	*        stops the run (see calc_q_par).
	*/
	void read_frames(const std::vector<std::string>& frame_paths,
		const Reference& ref, const Mesh& mesh, Background::Background& bkg,
		const bool strict)
	{
		const Zones::Zones& zones {bkg.get_zones()};
		const int nframes {static_cast<int>(frame_paths.size())};
		const int g {zones.get_ghost()};
		auto ne {Axisym::zoned(zones, nframes)};
		auto te {Axisym::zoned(zones, nframes)};
		auto ti {Axisym::zoned(zones, nframes)};
		auto vp {Axisym::zoned(zones, nframes)};
		auto upar {Axisym::zoned(zones, nframes)};
		std::vector<double> times {};

		// Parallel conductive heat flux (W/m^2) of each species, only if
		// the files have the face fluxes it comes from
		struct HeatFlux
		{
			std::string species;	// Group in the plasma files
			std::string name;		// Name in the background
			bool have {};
			Vectors::ZonedVector<BkgFPType> q_par {};
		};
		std::array<HeatFlux, 2> heat_fluxes {{{"spec0", "q_par_e"},
			{"spec1", "q_par_i"}}};

		for (int t {}; t < nframes; ++t)
		{
			const std::string& path {frame_paths[t]};
			const Hdf5::File file(path);
			std::cout << "Reading " << path << '\n';

			const int nzones {Hdf5::read_scalar<int>(file, "/NZones")};
			if (nzones != zones.get_nzones())
			{
				Errors::fatal_error(path + " has " + std::to_string(nzones)
					+ " zones, but the mesh has "
					+ std::to_string(zones.get_nzones()) + ".");
			}
			const auto nghosts {Hdf5::read_dataset<int>(file, "/Nghosts")};
			for (const int n : nghosts.data)
			{
				if (n != g)
				{
					Errors::fatal_error(path + " has " + std::to_string(n)
						+ " ghost cells on a side, but Flan expects "
						+ std::to_string(g) + ".");
				}
			}

			times.push_back(Hdf5::read_scalar<double>(file, "/tsimu")
				* ref.tau0);
			if (t == 0)
			{
				bkg.set_elec_mass_amu(Hdf5::read_scalar<double>(file,
					"/plasma_composition/spec0/mass"));
				bkg.set_ion_mass_amu(Hdf5::read_scalar<double>(file,
					"/plasma_composition/spec1/mass"));

			}
			for (auto& h : heat_fluxes)
			{
				const bool has {has_phi_face_fluxes(file, h.species)};
				if (t == 0)
				{
					h.have = has;
					if (has) h.q_par = Axisym::zoned(zones, nframes);
					else
					{
						std::cout << "  No " << h.species << " fluxes on the "
							"phi faces, so no heat flux (" << h.name
							<< ") is read.\n";
					}
				}
				else if (has != h.have)
				{
					Errors::fatal_error(path + (h.have ? " doesn't have"
						: " has") + " the " + h.species + " fluxes on the "
						"phi faces, but the first frame " + (h.have ? "does"
						: "doesn't") + ". Every frame needs them, or none.");
				}
				if (h.have)
				{
					const double mass {h.species == "spec0"
						? bkg.get_elec_mass_amu() : bkg.get_ion_mass_amu()};
					calc_q_par(file, ref, mesh, zones, h.species, mass, t,
						h.q_par, strict);
				}
			}

			for (int zone {}; zone < nzones; ++zone)
			{
				const std::string zg {zone_group(zone)};
				const std::vector<hsize_t> shape {
					static_cast<hsize_t>(zones.get_nx(zone) + 2 * g),
					static_cast<hsize_t>(zones.get_ny(zone) + 2 * g),
					static_cast<hsize_t>(zones.get_nz(zone) + 2 * g)};

				// The file's layout, ghost cells included, is the same as
				// this zone's block of a ZonedVector for frame t
				const std::size_t start {ne.calc_index(zone, t, -g, -g,
					-g)};
				auto load = [&](const std::string& name,
					Vectors::ZonedVector<BkgFPType>& v, const double scale)
				{
					const auto data {read_checked(file, zg + name, shape)};
					for (std::size_t n {}; n < data.size(); ++n)
						v.get_data()[start + n] = data[n] * scale;
				};
				load("/spec0/n", ne, ref.n0);
				load("/spec0/T", te, ref.T0);
				load("/spec1/T", ti, ref.T0);
				load("/PHI", vp, ref.phi0);

				const auto ni {read_checked(file, zg + "/spec1/n", shape)};
				const auto G {read_checked(file, zg + "/spec1/G", shape)};
				for (std::size_t n {}; n < ni.size(); ++n)
				{
					upar.get_data()[start + n] = ni[n] > 0.0
						? G[n] / ni[n] * ref.c0 : 0.0;
				}
			}

			if (t == 0)
			{
				std::cout << "Largest difference between SOLEDGE3X's ghost "
					"cells and Flan's, relative to the largest value, in "
					<< path << ":\n";
				Axisym::report_ghosts(zones, ne, "ne");
				Axisym::report_ghosts(zones, te, "Te");
				Axisym::report_ghosts(zones, ti, "Ti");
				Axisym::report_ghosts(zones, vp, "vp");
				Axisym::report_ghosts(zones, upar, "u_par");
			}
		}

		for (int t {1}; t < nframes; ++t)
		{
			if (!(times[t] > times[t - 1]))
			{
				Errors::fatal_error("The SOLEDGE3X frames' times don't "
					"increase: " + frame_paths[t] + " is at "
					+ std::to_string(times[t]) + " s, but the frame before "
					"it is at " + std::to_string(times[t - 1]) + " s.");
			}
		}

		for (auto* v : {&ne, &te, &ti, &vp, &upar}) zones.fill_ghosts(*v);

		// Ti goes in first, since the flows need it
		bkg.move_into_times(times);
		bkg.move_into_ne(ne);
		bkg.move_into_te(te);
		bkg.move_into_ti(ti);
		for (auto& h : heat_fluxes)
		{
			if (!h.have) continue;
			zones.fill_ghosts(h.q_par);
			if (h.species == "spec0") bkg.move_into_q_par_e(h.q_par);
			else bkg.move_into_q_par_i(h.q_par);
		}
		Axisym::calc_efield_and_flows(bkg, vp, upar);
		bkg.move_into_vp(vp);
	}

	/**
	* Reads a SOLEDGE3X simulation into a Background.
	*
	* The files are in s3x_dir: the mesh (s3x_mesh_file), the reference
	* values (s3x_ref_file) and the plasma frames plasma_%05d.h5 from
	* s3x_frame_start to s3x_frame_end, of which there must be at least
	* two. s3x_dir must be set. The reference values are printed.
	*
	* @param opts Simulation options.
	* @return The background.
	*/
	Background::Background read_s3xe(const Options::Options& opts)
	{
		const std::string dir {opts.s3x_dir()};
		if (dir == "undefined")
		{
			Errors::fatal_error("s3x_dir must be set to the directory with "
				"the SOLEDGE3X files.");
		}
		const int start {opts.s3x_frame_start()};
		const int end {opts.s3x_frame_end()};
		if (end - start + 1 < 2)
		{
			Errors::fatal_error("A SOLEDGE3X background needs at least two "
				"frames, but s3x_frame_start = " + std::to_string(start)
				+ " and s3x_frame_end = " + std::to_string(end) + ".");
		}
		std::vector<std::string> frame_paths {};
		for (int frame {start}; frame <= end; ++frame)
		{
			char name[32];
			std::snprintf(name, sizeof(name), "plasma_%05d.h5", frame);
			frame_paths.push_back(dir + "/" + name);
		}

		const Reference ref {read_reference(dir + "/" + opts.s3x_ref_file())};
		print_reference(ref);

		const int nphi {read_nphi(frame_paths[0])};
		const Mesh mesh {read_mesh(dir + "/" + opts.s3x_mesh_file(), nphi)};
		check_reference(ref, mesh);

		Background::Background bkg {};
		auto zones {mesh.zones};
		bkg.move_into_zones(zones);
		{
			const Axisym::Mesh geom {to_axisym(mesh)};
			Axisym::calc_geometry(geom, bkg);
			Axisym::calc_bfield(geom, static_cast<int>(frame_paths.size()),
				bkg);
		}
		load_chi(mesh, bkg);
		bkg.set_walls(mesh.walls);

		// The ion temperature is needed for the grad B drift, so the
		// plasma (which includes calculating the flows) comes last
		read_frames(frame_paths, ref, mesh, bkg,
			opts.imp_thermal_force_int() == 1);

		return bkg;
	}
}
