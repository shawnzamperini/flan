/*
* Writes small SOLEDGE3X-style mesh and plasma files for tests, so they
* don't depend on real (multi-GB) SOLEDGE3X output.
*
* The test mesh has three zones on a circular cross-section, in normalized
* units R = R0 + a0 * r cos(theta) and Z = a0 * r sin(theta):
*
* - Zones 1 and 2 are a ring like the core: r from r_in to r_mid, theta
*   from 0 to pi (zone 1) and pi to 2 pi (zone 2), connected to each
*   other at both theta ends. Their psi- faces are the core.
* - Zone 3 is outside zone 1 in psi: r from r_mid to r_out, theta from 0
*   to pi, with targets (-1) at both theta ends and the outer boundary
*   (-3) at psi+. Zone 2's psi+ is also the outer boundary.
*
* Br = Bz = 0 and Bphi = R0 / R, a vacuum toroidal field. chi is 1 in the
* last psi row of zone 3 and 0 elsewhere. The wall outline (/walls) is a
* circle along the inner edge of that row, so it runs along the faces
* between zone 3's plasma and wall cells.
*/

#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "hdf5.h"
#include "zones.h"

namespace S3XTest
{
	constexpr double pi {3.141592653589793};

	// Normalization and layout of the test mesh
	struct TestMesh
	{
		double R0 {1.7};
		double a0 {0.5};
		double B0 {2.0};
		double r_in {0.2};
		double r_mid {0.5};
		double r_out {0.7};
		int npsi_core {12};		// Cells in psi, zones 1 and 2
		int npsi_sol {8};		// Cells in psi, zone 3
		int ntheta {40};		// Cells in theta, every zone
		int nphi {6};			// Toroidal planes, in the plasma file

		// Things the tests can break
		int nzones_written {3};		// Zones written to NZones
		int neighbors_cols {3};		// Columns of the neighbors table
		int bad_code {0};			// If not 0, zone 3's psi+ code
		int mesh_nphi {1};			// Nphi written in the mesh file
		bool skip_zone3 {false};	// Leave out zone 3's group
		int corner_extra {0};		// Extra rows in zone 1's Rcorners

		// Wall options
		bool write_walls {true};	// Write the /walls outline
		int wall_points {400};		// Points around the wall circle
		bool chi_zone3_first_row {false};	// Also chi = 1 in zone 3's
											// first psi row
	};

	// Normalized r of the wall outline: the inner edge of zone 3's last psi
	// row, where chi = 1
	inline double wall_r(const TestMesh& m)
	{
		return m.r_mid + (m.r_out - m.r_mid) * (m.npsi_sol - 1) / m.npsi_sol;
	}

	// Inner and outer r and the theta range of a zone (0-based)
	inline void zone_extent(const TestMesh& m, int zone, double& r0,
		double& r1, double& th0, double& th1, int& npsi)
	{
		r0 = zone == 2 ? m.r_mid : m.r_in;
		r1 = zone == 2 ? m.r_out : m.r_mid;
		th0 = zone == 1 ? pi : 0.0;
		th1 = zone == 1 ? 2.0 * pi : pi;
		npsi = zone == 2 ? m.npsi_sol : m.npsi_core;
	}

	// Normalized r and theta of corner (i, j) of a zone
	inline double corner_r(const TestMesh& m, int zone, int i)
	{
		double r0, r1, th0, th1;
		int npsi;
		zone_extent(m, zone, r0, r1, th0, th1, npsi);
		return r0 + (r1 - r0) * i / npsi;
	}

	inline double corner_theta(const TestMesh& m, int zone, int j)
	{
		double r0, r1, th0, th1;
		int npsi;
		zone_extent(m, zone, r0, r1, th0, th1, npsi);
		return th0 + (th1 - th0) * j / m.ntheta;
	}

	// R (m) at the bilinear center of cell (i, j), the same way Flan
	// calculates it from the corners
	inline double cell_R(const TestMesh& m, int zone, int i, int j)
	{
		double sum {};
		for (int di {}; di < 2; ++di)
		for (int dj {}; dj < 2; ++dj)
		{
			sum += m.R0 + m.a0 * corner_r(m, zone, i + di)
				* std::cos(corner_theta(m, zone, j + dj));
		}
		return 0.25 * sum;
	}

	// Z (m) at the bilinear center of cell (i, j)
	inline double cell_Z(const TestMesh& m, int zone, int i, int j)
	{
		double sum {};
		for (int di {}; di < 2; ++di)
		for (int dj {}; dj < 2; ++dj)
		{
			sum += m.a0 * corner_r(m, zone, i + di)
				* std::sin(corner_theta(m, zone, j + dj));
		}
		return 0.25 * sum;
	}

	// Normalized plasma of the test frames, as functions of a cell's
	// center (R, Z) in meters and the frame number t
	inline double test_ne(double R, double Z, int t)
		{return 2.0 + Z + 0.5 * t + 0.1 * R;}
	inline double test_te(double R, double, int t)
		{return 0.5 + 0.1 * R + 0.05 * t;}
	inline double test_ti(double R, double, int t)
		{return 0.8 + 0.1 * R + 0.05 * t;}
	// PHI is linear in Z, so E = -0.5 phi0 Z_hat
	inline double test_phi(double, double Z, int) {return 0.5 * Z;}
	// G / n = 0.3, so u_par = 0.3 c0
	inline double test_G(double R, double Z, int t)
		{return 0.3 * test_ne(R, Z, t);}
	// Ion parallel conductive heat flux (normalized by n0 c0 T0), positive
	// along B, written into the phi-face fluxes when a frame has them
	inline double test_q_par(double, double Z, int t)
		{return 0.05 * Z - 0.02 + 0.01 * t;}
	// Electron G / n = 0.2, and the electron heat flux like test_q_par
	inline double test_G_e(double R, double Z, int t)
		{return 0.2 * test_ne(R, Z, t);}
	inline double test_q_par_e(double, double Z, int t)
		{return -0.03 * Z + 0.5 + 0.02 * t;}

	// Poloidal area of a cell of the test mesh, normalized by a0^2, from
	// the shoelace formula on its four corners
	inline double cell_area(const TestMesh& m, int zone, int i, int j)
	{
		const std::array<std::array<int, 2>, 4> c {{{i, j}, {i + 1, j},
			{i + 1, j + 1}, {i, j + 1}}};
		double twice {};
		for (int p {}; p < 4; ++p)
		{
			const auto& a {c[p]};
			const auto& b {c[(p + 1) % 4]};
			const double ra {corner_r(m, zone, a[0])};
			const double ta {corner_theta(m, zone, a[1])};
			const double rb {corner_r(m, zone, b[0])};
			const double tb {corner_theta(m, zone, b[1])};
			twice += ra * std::cos(ta) * rb * std::sin(tb)
				- rb * std::cos(tb) * ra * std::sin(ta);
		}
		return 0.5 * std::abs(twice);
	}

	inline void write(hid_t file, const std::string& name,
		const std::vector<hsize_t>& shape, const std::vector<double>& data)
	{
		const hid_t lcpl {H5Pcreate(H5P_LINK_CREATE)};
		H5Pset_create_intermediate_group(lcpl, 1);
		const hid_t space {H5Screate_simple(static_cast<int>(shape.size()),
			shape.data(), nullptr)};
		const hid_t dset {H5Dcreate2(file, name.c_str(), H5T_IEEE_F64LE,
			space, lcpl, H5P_DEFAULT, H5P_DEFAULT)};
		H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			data.data());
		H5Dclose(dset);
		H5Sclose(space);
		H5Pclose(lcpl);
	}

	// Integers are written as 64-bit, like SOLEDGE3X does
	inline void write(hid_t file, const std::string& name,
		const std::vector<hsize_t>& shape, const std::vector<int64_t>& data)
	{
		const hid_t lcpl {H5Pcreate(H5P_LINK_CREATE)};
		H5Pset_create_intermediate_group(lcpl, 1);
		const hid_t space {H5Screate_simple(static_cast<int>(shape.size()),
			shape.data(), nullptr)};
		const hid_t dset {H5Dcreate2(file, name.c_str(), H5T_STD_I64LE,
			space, lcpl, H5P_DEFAULT, H5P_DEFAULT)};
		H5Dwrite(dset, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			data.data());
		H5Dclose(dset);
		H5Sclose(space);
		H5Pclose(lcpl);
	}

	inline void write_mesh(const std::string& path, const TestMesh& m)
	{
		const hid_t file {H5Fcreate(path.c_str(), H5F_ACC_TRUNC,
			H5P_DEFAULT, H5P_DEFAULT)};
		write(file, "/R0", {1}, std::vector<double> {m.R0});
		write(file, "/a0", {1}, std::vector<double> {m.a0});
		write(file, "/B0", {1}, std::vector<double> {m.B0});
		write(file, "/NZones", {1},
			std::vector<int64_t> {m.nzones_written});

		// Rows psi-, psi+, theta-, theta+, phi-, phi+; columns zones 1-3
		std::vector<int64_t> nb {
			-2, -2, 1,
			3, -3, m.bad_code ? m.bad_code : -3,
			2, 1, -1,
			2, 1, -1,
			1, 2, 3,
			1, 2, 3};
		if (m.neighbors_cols != 3)
		{
			std::vector<int64_t> resized {};
			for (int row {}; row < 6; ++row)
			for (int col {}; col < m.neighbors_cols; ++col)
				resized.push_back(col < 3 ? nb[row * 3 + col] : -3);
			nb = resized;
		}
		write(file, "/neighbors",
			{6, static_cast<hsize_t>(m.neighbors_cols)}, nb);

		for (int zone {}; zone < 3; ++zone)
		{
			if (zone == 2 && m.skip_zone3) continue;
			double r0, r1, th0, th1;
			int npsi;
			zone_extent(m, zone, r0, r1, th0, th1, npsi);
			const int nt {m.ntheta};
			const std::string g {"/zone" + std::to_string(zone + 1)};

			write(file, g + "/Npsi", {1}, std::vector<int64_t> {npsi});
			write(file, g + "/Ntheta", {1}, std::vector<int64_t> {nt});
			write(file, g + "/Nphi", {1},
				std::vector<int64_t> {m.mesh_nphi});
			write(file, g + "/phicorners", {2},
				std::vector<double> {0.0, 2.0 * pi});

			const int extra {zone == 0 ? m.corner_extra : 0};
			std::vector<double> Rc {}, Zc {};
			for (int i {}; i <= npsi + extra; ++i)
			for (int j {}; j <= nt; ++j)
			{
				const double r {corner_r(m, zone, i)};
				const double th {corner_theta(m, zone, j)};
				Rc.push_back(r * std::cos(th));
				Zc.push_back(r * std::sin(th));
			}
			const hsize_t np {static_cast<hsize_t>(npsi)};
			const hsize_t nth {static_cast<hsize_t>(nt)};
			write(file, g + "/Rcorners",
				{np + 1 + static_cast<hsize_t>(extra), nth + 1, 1}, Rc);
			write(file, g + "/Zcorners",
				{np + 1 + static_cast<hsize_t>(extra), nth + 1, 1}, Zc);

			std::vector<double> zero(npsi * nt, 0.0), bphi {}, chi {};
			for (int i {}; i < npsi; ++i)
			for (int j {}; j < nt; ++j)
			{
				bphi.push_back(m.R0 / cell_R(m, zone, i, j));
				const bool wall {zone == 2 && (i == npsi - 1
					|| (m.chi_zone3_first_row && i == 0))};
				chi.push_back(wall ? 1.0 : 0.0);
			}
			write(file, g + "/Br", {np, nth, 1}, zero);
			write(file, g + "/Bz", {np, nth, 1}, zero);
			write(file, g + "/Bphi", {np, nth, 1}, bphi);
			write(file, g + "/chi", {np, nth, 1}, chi);
		}

		// Wall outline in meters (not normalized), clockwise in (R, Z) and
		// closed by repeating the first point, like SOLEDGE3X's
		if (m.write_walls)
		{
			std::vector<double> wR {}, wZ {};
			for (int p {}; p <= m.wall_points; ++p)
			{
				const double th {-2.0 * pi * p / m.wall_points};
				wR.push_back(m.R0 + m.a0 * wall_r(m) * std::cos(th));
				wZ.push_back(m.a0 * wall_r(m) * std::sin(th));
			}
			const hsize_t n {static_cast<hsize_t>(wR.size())};
			write(file, "/walls/Nwalls", {1}, std::vector<double> {1.0});
			write(file, "/walls/wall1/R", {n}, wR);
			write(file, "/walls/wall1/Z", {n}, wZ);
			write(file, "/walls/wall1/type", {1}, std::vector<double> {1.0});
		}
		H5Fclose(file);
	}

	// Reference file with the same values as the DIII-D case, except R0, a0
	// and B0 from the test mesh (or B0_ref if given), and T0 = T0_ref
	inline void write_reference(const std::string& path, const TestMesh& m,
		double B0_ref = 0.0, double T0_ref = 100.0)
	{
		const hid_t file {H5Fcreate(path.c_str(), H5F_ACC_TRUNC,
			H5P_DEFAULT, H5P_DEFAULT)};
		write(file, "/n0", {1}, std::vector<double> {1e19});
		write(file, "/T0", {1}, std::vector<double> {T0_ref});
		write(file, "/c0", {1}, std::vector<double> {98227.0});
		write(file, "/tau0", {1}, std::vector<double> {5.23951e-9});
		write(file, "/phi0", {1}, std::vector<double> {100.0});
		write(file, "/R0", {1}, std::vector<double> {m.R0});
		write(file, "/a0", {1}, std::vector<double> {m.a0});
		write(file, "/B0", {1},
			std::vector<double> {B0_ref > 0.0 ? B0_ref : m.B0});
		H5Fclose(file);
	}

	// Things a test can break in a plasma file
	struct TestFrame
	{
		int nzones_written {3};
		int nghosts {2};
		double tsimu {};
		int extra_theta {0};	// Extra theta cells in zone 2's datasets

		// If true, also write the electron and ion fluxes through the phi
		// faces, with fluxn = flux_C G b_phi A_pol / a0^2 and fluxE =
		// (5/2 T + 1/2 m u^2) fluxn + flux_C q_par b_phi A_pol / a0^2
		// (b_phi = 1 here)
		bool face_fluxes {false};
		double flux_C {-0.99 * 2.0 * pi * 98227.0 / 0.5};
	};

	// Plasma file for frame t of the test mesh, whose zones (from
	// S3X::read_mesh) give the ghost cells. Ghost cells get the value of
	// the cells Flan fills them from, so they match fill_ghosts exactly.
	inline void write_plasma(const std::string& path, const TestMesh& m,
		const Zones::Zones& zones, int t, const TestFrame& frame)
	{
		const hid_t file {H5Fcreate(path.c_str(), H5F_ACC_TRUNC,
			H5P_DEFAULT, H5P_DEFAULT)};
		const int g {2};
		write(file, "/NZones", {1},
			std::vector<int64_t> {frame.nzones_written});
		write(file, "/Nghosts", {6}, std::vector<int64_t>(6, frame.nghosts));
		write(file, "/tsimu", {1}, std::vector<double> {frame.tsimu});
		write(file, "/plasma_composition/spec0/mass", {1},
			std::vector<double> {0.00054858});
		write(file, "/plasma_composition/spec1/mass", {1},
			std::vector<double> {2.014});

		using Field = double (*)(double, double, int);
		const std::vector<std::pair<std::string, Field>> fields {
			{"/spec0/n", test_ne}, {"/spec0/T", test_te},
			{"/spec1/T", test_ti}, {"/spec1/n", test_ne},
			{"/spec1/G", test_G}, {"/spec0/G", test_G_e},
			{"/PHI", test_phi}};

		for (int zone {}; zone < 3; ++zone)
		{
			const int nx {zones.get_nx(zone)};
			const int ny {zones.get_ny(zone)
				+ (zone == 1 ? frame.extra_theta : 0)};
			const int nz {zones.get_nz(zone)};
			const std::string zg {"/zone" + std::to_string(zone + 1)};
			write(file, zg + "/Nphi", {1}, std::vector<int64_t> {nz});

			for (const auto& [name, f] : fields)
			{
				std::vector<double> data {};
				for (int i {-g}; i < nx + g; ++i)
				for (int j {-g}; j < ny + g; ++j)
				for (int k {-g}; k < nz + g; ++k)
				{
					// Average over the cells this one is filled from
					const auto sources {zones.resolve_cell({zone, i,
						std::min(j, zones.get_ny(zone) + g - 1), k})};
					double sum {};
					for (const auto& c : sources)
					{
						sum += f(cell_R(m, c.zone, c.i, c.j),
							cell_Z(m, c.zone, c.i, c.j), t);
					}
					data.push_back(sum / sources.size());
				}
				write(file, zg + name, {static_cast<hsize_t>(nx + 2 * g),
					static_cast<hsize_t>(ny + 2 * g),
					static_cast<hsize_t>(nz + 2 * g)}, data);
			}

			// Nothing depends on phi, so the face values are the cell's
			for (const bool elec : {false, true})
			{
				if (!frame.face_fluxes) break;
				const Field Gf {elec ? test_G_e : test_G};
				const Field Tf {elec ? test_te : test_ti};
				const Field qf {elec ? test_q_par_e : test_q_par};
				const double mass {elec ? 0.00054858 : 2.014};
				std::vector<double> fluxn {}, fluxE {};
				for (int i {}; i < nx; ++i)
				for (int j {}; j < ny; ++j)
				for (int k {}; k <= nz; ++k)
				{
					const double R {cell_R(m, zone, i, j)};
					const double Z {cell_Z(m, zone, i, j)};
					const double n {test_ne(R, Z, t)};
					const double G {Gf(R, Z, t)};
					const double u {G / n};
					const double h {2.5 * Tf(R, Z, t) + 0.5 * mass * u * u};
					const double geom {frame.flux_C
						* cell_area(m, zone, i, j)};
					fluxn.push_back(geom * G);
					fluxE.push_back(h * geom * G + geom * qf(R, Z, t));
				}
				const std::vector<hsize_t> shape {static_cast<hsize_t>(nx),
					static_cast<hsize_t>(ny), static_cast<hsize_t>(nz + 1)};
				const std::string sp {elec ? "/spec0" : "/spec1"};
				write(file, zg + sp + "/fluxn/phi", shape, fluxn);
				write(file, zg + sp + "/fluxE/phi", shape, fluxE);
			}
		}
		H5Fclose(file);
	}

	// Plasma file with only each zone's Nphi, for reading the number of
	// toroidal planes. zone3_nphi lets a test make one zone different.
	inline void write_nphi_only(const std::string& path, const TestMesh& m,
		int zone3_nphi)
	{
		const hid_t file {H5Fcreate(path.c_str(), H5F_ACC_TRUNC,
			H5P_DEFAULT, H5P_DEFAULT)};
		write(file, "/NZones", {1}, std::vector<int64_t> {3});
		for (int zone {}; zone < 3; ++zone)
		{
			write(file, "/zone" + std::to_string(zone + 1) + "/Nphi", {1},
				std::vector<int64_t> {zone == 2 ? zone3_nphi : m.nphi});
		}
		H5Fclose(file);
	}

	// Writes a whole test case to dir, with the names read_s3xe expects:
	// mesh_S3X.h5, refParam_raptorX.h5, and plasma_00000.h5 and
	// plasma_00001.h5 at tsimu 1e7 and 1e7 + 1000. zones must be the test
	// mesh's zones (from S3X::read_mesh), which the plasma ghost cells come
	// from, so it is passed in instead of read here.
	template <typename ReadMesh>
	void write_case(const std::string& dir, const TestMesh& m,
		ReadMesh read_mesh)
	{
		write_mesh(dir + "/mesh_S3X.h5", m);
		write_reference(dir + "/refParam_raptorX.h5", m);
		const auto zones {read_mesh(dir + "/mesh_S3X.h5")};
		for (int t {}; t < 2; ++t)
		{
			TestFrame frame {};
			frame.tsimu = 1.0e7 + 1000.0 * t;
			char name[32];
			std::snprintf(name, sizeof(name), "/plasma_%05d.h5", t);
			write_plasma(dir + name, m, zones, t, frame);
		}
	}
}
