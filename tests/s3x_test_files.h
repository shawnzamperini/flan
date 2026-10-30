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
* last psi row of zone 3 and 0 elsewhere.
*/

#pragma once

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "hdf5.h"

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
	};

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
				chi.push_back(zone == 2 && i == npsi - 1 ? 1.0 : 0.0);
			}
			write(file, g + "/Br", {np, nth, 1}, zero);
			write(file, g + "/Bz", {np, nth, 1}, zero);
			write(file, g + "/Bphi", {np, nth, 1}, bphi);
			write(file, g + "/chi", {np, nth, 1}, chi);
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

	// Plasma file with only each zone's Nphi, for reading the number of
	// toroidal planes. zone3_nphi lets a test make one zone different.
	inline void write_nphi_only(const std::string& path, const TestMesh& m,
		int zone3_nphi)
	{
		const hid_t file {H5Fcreate(path.c_str(), H5F_ACC_TRUNC,
			H5P_DEFAULT, H5P_DEFAULT)};
		for (int zone {}; zone < 3; ++zone)
		{
			write(file, "/zone" + std::to_string(zone + 1) + "/Nphi", {1},
				std::vector<int64_t> {zone == 2 ? zone3_nphi : m.nphi});
		}
		H5Fclose(file);
	}
}
