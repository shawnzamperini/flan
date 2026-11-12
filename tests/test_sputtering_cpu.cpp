/*
* Tests for sputtering:
*
* - Reading and interpolating sputtering yield tables
*   (Sputtering::YieldTable), on small tables written by the test itself
*   and on the table in data/sputtering/.
* - Finding the wall surfaces (Sputtering::Surface), on a small background
*   built by hand with a wall at 45 degrees to the grid, where every answer
*   is known exactly, and on the SOLEDGE3X test mesh of s3x_test_files.h,
*   whose wall is a circle.
*/

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "background.h"
#include "flan_types.h"
#include "hdf5.h"
#include "hdf5_utils.h"
#include "read_s3xe.h"
#include "s3x_test_files.h"
#include "sputtering.h"
#include "vectors.h"
#include "zones.h"


namespace
{
	// Writes a dataset of doubles with the given shape (empty for a
	// scalar), creating groups along the way
	void write_doubles(hid_t file, const std::string& name,
		const std::vector<hsize_t>& shape, const std::vector<double>& data)
	{
		const hid_t lcpl {H5Pcreate(H5P_LINK_CREATE)};
		H5Pset_create_intermediate_group(lcpl, 1);
		const hid_t space {shape.empty() ? H5Screate(H5S_SCALAR)
			: H5Screate_simple(static_cast<int>(shape.size()), shape.data(),
			nullptr)};
		const hid_t dset {H5Dcreate2(file, name.c_str(), H5T_NATIVE_DOUBLE,
			space, lcpl, H5P_DEFAULT, H5P_DEFAULT)};
		H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			data.data());
		H5Dclose(dset);
		H5Sclose(space);
		H5Pclose(lcpl);
	}

	// The test table's yield above its lowest energy. Linear in log10(E)
	// and angle, so bilinear interpolation gives it exactly.
	double linear_yield(double energy, double angle)
	{
		return 0.01 + 0.02 * std::log10(energy) + 0.001 * angle;
	}

	// What goes in a test table. The defaults are a good table: 4 energies
	// and 4 angles, with the lowest energy below threshold (Y = 0) and
	// linear_yield above it.
	struct Table
	{
		std::vector<double> energy {10.0, 100.0, 1000.0, 10000.0};
		std::vector<double> angle {0.0, 30.0, 60.0, 89.0};
		std::vector<double> yield {};
		std::vector<hsize_t> yield_shape {};	// Default (nE, nA)
		bool target_Es {true};					// Write /target/Es

		Table()
		{
			for (const double e : energy)
			{
				for (const double a : angle)
				{
					yield.push_back(e == 10.0 ? 0.0 : linear_yield(e, a));
				}
			}
		}
	};

	// Test tables are written here and removed when the test ends
	class SputteringTable : public testing::Test
	{
	protected:
		std::string m_path {(std::filesystem::temp_directory_path()
			/ ("flan_test_sputtering_" + std::to_string(getpid())
			+ ".h5")).string()};

		void write(const Table& t)
		{
			const hid_t file {H5Fcreate(m_path.c_str(), H5F_ACC_TRUNC,
				H5P_DEFAULT, H5P_DEFAULT)};
			write_doubles(file, "/energy", {t.energy.size()}, t.energy);
			write_doubles(file, "/angle", {t.angle.size()}, t.angle);
			write_doubles(file, "/yield", t.yield_shape.empty()
				? std::vector<hsize_t> {t.energy.size(), t.angle.size()}
				: t.yield_shape, t.yield);
			write_doubles(file, "/projectile/Z", {}, {1.0});
			write_doubles(file, "/projectile/mass_amu", {}, {2.014});
			write_doubles(file, "/target/Z", {}, {74.0});
			write_doubles(file, "/target/mass_amu", {}, {183.84});
			if (t.target_Es) write_doubles(file, "/target/Es", {}, {8.79});
			H5Fclose(file);
		}

		void TearDown() override {std::filesystem::remove(m_path);}
	};
}

// The table's values come back exactly at the grid points, and the
// species and energy range are read
TEST_F(SputteringTable, Read)
{
	write(Table {});
	const Sputtering::YieldTable table(m_path);

	EXPECT_EQ(table.path(), m_path);
	EXPECT_DOUBLE_EQ(table.energy_min(), 10.0);
	EXPECT_DOUBLE_EQ(table.energy_max(), 10000.0);
	EXPECT_DOUBLE_EQ(table.projectile_Z(), 1.0);
	EXPECT_DOUBLE_EQ(table.projectile_mass_amu(), 2.014);
	EXPECT_DOUBLE_EQ(table.target_Z(), 74.0);
	EXPECT_DOUBLE_EQ(table.target_mass_amu(), 183.84);
	EXPECT_DOUBLE_EQ(table.target_Es(), 8.79);

	const Table t {};
	for (std::size_t i {}; i < t.energy.size(); ++i)
	{
		for (std::size_t j {}; j < t.angle.size(); ++j)
		{
			EXPECT_NEAR(table.get_yield(t.energy[i], t.angle[j]),
				t.yield[i * t.angle.size() + j], 1e-15);
		}
	}
}

// Between grid points the yield is linear in log10(E) and angle
TEST_F(SputteringTable, Interpolates)
{
	write(Table {});
	const Sputtering::YieldTable table(m_path);

	for (const double e : {150.0, 316.2, 999.0, 2500.0, 7777.0})
	{
		for (const double a : {0.0, 12.5, 45.0, 61.0, 88.0})
		{
			EXPECT_NEAR(table.get_yield(e, a), linear_yield(e, a), 1e-14)
				<< "E = " << e << ", angle = " << a;
		}
	}

	// Between the zero row and the next one, the yield goes linearly in
	// log10(E) from 0, e.g., half way at the geometric mean
	const double e_mid {std::sqrt(10.0 * 100.0)};
	EXPECT_NEAR(table.get_yield(e_mid, 30.0),
		0.5 * linear_yield(100.0, 30.0), 1e-14);
}

// Below the table the yield is 0, above it the highest energy is used,
// and angles are clamped to the table's
TEST_F(SputteringTable, OutsideTable)
{
	write(Table {});
	const Sputtering::YieldTable table(m_path);

	EXPECT_EQ(table.get_yield(9.99, 30.0), 0.0);
	EXPECT_EQ(table.get_yield(0.0, 30.0), 0.0);
	EXPECT_EQ(table.get_yield(-5.0, 30.0), 0.0);
	EXPECT_EQ(table.get_yield(NAN, 30.0), 0.0);

	EXPECT_NEAR(table.get_yield(1e6, 30.0), linear_yield(10000.0, 30.0),
		1e-14);
	EXPECT_NEAR(table.get_yield(1000.0, 90.0), linear_yield(1000.0, 89.0),
		1e-14);
	EXPECT_NEAR(table.get_yield(1000.0, -10.0), linear_yield(1000.0, 0.0),
		1e-14);
}

// Bad tables stop the run with a message saying what's wrong
TEST_F(SputteringTable, BadTables)
{
	EXPECT_EXIT(Sputtering::YieldTable("/no/such/table.h5"),
		testing::ExitedWithCode(EXIT_FAILURE), "table.h5");

	Table missing {};
	missing.target_Es = false;
	write(missing);
	EXPECT_EXIT(Sputtering::YieldTable {m_path},
		testing::ExitedWithCode(EXIT_FAILURE), "/target/Es");

	Table shape {};
	shape.yield_shape = {4, 2, 2};
	write(shape);
	EXPECT_EXIT(Sputtering::YieldTable {m_path},
		testing::ExitedWithCode(EXIT_FAILURE), "wrong shape");

	Table one {};
	one.energy = {100.0};
	one.yield = {0.1, 0.1, 0.1, 0.1};
	write(one);
	EXPECT_EXIT(Sputtering::YieldTable {m_path},
		testing::ExitedWithCode(EXIT_FAILURE), "At least 2");

	Table order {};
	order.energy = {10.0, 1000.0, 100.0, 10000.0};
	write(order);
	EXPECT_EXIT(Sputtering::YieldTable {m_path},
		testing::ExitedWithCode(EXIT_FAILURE), "not strictly increasing");

	Table zero {};
	zero.energy[0] = 0.0;
	write(zero);
	EXPECT_EXIT(Sputtering::YieldTable {m_path},
		testing::ExitedWithCode(EXIT_FAILURE), "energies <= 0");

	Table angle {};
	angle.angle.back() = 95.0;
	write(angle);
	EXPECT_EXIT(Sputtering::YieldTable {m_path},
		testing::ExitedWithCode(EXIT_FAILURE), "outside 0 to 90");

	Table negative {};
	negative.yield[5] = -0.1;
	write(negative);
	EXPECT_EXIT(Sputtering::YieldTable {m_path},
		testing::ExitedWithCode(EXIT_FAILURE), "negative or non-finite");
}

// The D on W table in the repository reads, and gives its own values back
// at its grid points
TEST(SputteringData, DOnW)
{
	const std::string path {std::string {FLAN_SOURCE_DIR}
		+ "/data/sputtering/D_on_W.h5"};
	if (!std::filesystem::exists(path))
	{
		GTEST_SKIP() << path << " doesn't exist yet.";
	}

	const Sputtering::YieldTable table(path);
	EXPECT_DOUBLE_EQ(table.projectile_Z(), 1.0);
	EXPECT_DOUBLE_EQ(table.target_Z(), 74.0);

	const Hdf5::File file(path);
	const auto energy {Hdf5::read_dataset<double>(file, "/energy").data};
	const auto angle {Hdf5::read_dataset<double>(file, "/angle").data};
	const auto yield {Hdf5::read_dataset<double>(file, "/yield").data};
	for (std::size_t i {}; i < energy.size(); ++i)
	{
		for (std::size_t j {}; j < angle.size(); ++j)
		{
			EXPECT_NEAR(table.get_yield(energy[i], angle[j]),
				yield[i * angle.size() + j],
				1e-12 * yield[i * angle.size() + j])
				<< "E = " << energy[i] << ", angle = " << angle[j];
		}
	}
}


namespace
{
	// The background for the 45 degree wall tests: one zone of 10 x 10
	// cells, each 1 cm in R and Z, starting at R = 1 m and Z = 0, and 4
	// toroidal cells of 0.1 rad. x is R and y is Z, both as cell indices.
	//
	// Cells with i + j >= 10 are wall (chi = 1), so the plasma is the lower
	// left triangle in (R, Z) and the wall the upper right, with a
	// staircase between them. Cell (2, 2) is wall too, alone in the
	// plasma, to stand for wall the outline doesn't include. The outline
	// is the straight line x + y = 10.5 (in cell units), which goes
	// through the middle of every face of the staircase.
	constexpr int diag_n {10};
	constexpr int diag_nz {4};
	constexpr double diag_dphi {0.1};

	double diag_R(double x) {return 1.0 + 0.01 * x;}
	double diag_Z(double y) {return 0.01 * y;}

	Background::Background diagonal_background(bool walls = true)
	{
		const int n {diag_n};
		const int nz {diag_nz};
		std::vector<double> x(n), gx(n + 1), z(nz), gz(nz + 1);
		for (int i {}; i <= n; ++i) gx[i] = i;
		for (int i {}; i < n; ++i) x[i] = i + 0.5;
		for (int k {}; k <= nz; ++k) gz[k] = k * diag_dphi;
		for (int k {}; k < nz; ++k) z[k] = (k + 0.5) * diag_dphi;

		Zones::Zones zones(2);
		zones.add_zone(x, x, z, gx, gx, gz);
		Background::Background bkg {};
		std::vector<BkgFPType> times {0.0, 1.0};
		bkg.move_into_times(times);
		bkg.move_into_zones(zones);

		// Cell centers and corners, in Cartesian coordinates
		Vectors::ZonedVector<BkgFPType> X(1, {n}, {n}, {nz}, 2);
		Vectors::ZonedVector<BkgFPType> Y(1, {n}, {n}, {nz}, 2);
		Vectors::ZonedVector<BkgFPType> Z(1, {n}, {n}, {nz}, 2);
		Vectors::ZonedVector<BkgFPType> chi(1, {n}, {n}, {nz}, 2);
		for (int i {}; i < n; ++i)
		for (int j {}; j < n; ++j)
		for (int k {}; k < nz; ++k)
		{
			X(0, i, j, k) = diag_R(i + 0.5) * std::cos(z[k]);
			Y(0, i, j, k) = diag_R(i + 0.5) * std::sin(z[k]);
			Z(0, i, j, k) = diag_Z(j + 0.5);
			chi(0, i, j, k) = (i + j >= n || (i == 2 && j == 2)) ? 1.0
				: 0.0;
		}
		bkg.get_zones().fill_ghosts(chi);

		Vectors::ZonedVector<BkgFPType> gX(1, {n + 1}, {n + 1}, {nz + 1}, 0);
		Vectors::ZonedVector<BkgFPType> gY(1, {n + 1}, {n + 1}, {nz + 1}, 0);
		Vectors::ZonedVector<BkgFPType> gZ(1, {n + 1}, {n + 1}, {nz + 1}, 0);
		for (int i {}; i <= n; ++i)
		for (int j {}; j <= n; ++j)
		for (int k {}; k <= nz; ++k)
		{
			gX(0, i, j, k) = diag_R(i) * std::cos(gz[k]);
			gY(0, i, j, k) = diag_R(i) * std::sin(gz[k]);
			gZ(0, i, j, k) = diag_Z(j);
		}

		bkg.move_into_X(X);
		bkg.move_into_Y(Y);
		bkg.move_into_Z(Z);
		bkg.move_into_grid_X(gX);
		bkg.move_into_grid_Y(gY);
		bkg.move_into_grid_Z(gZ);
		bkg.move_into_chi(chi);

		// The outline x + y = 10.5, from the R axis to the Z axis
		if (walls)
		{
			Background::WallPolygon wall {};
			wall.R = {diag_R(10.5), diag_R(0.0)};
			wall.Z = {diag_Z(0.0), diag_Z(10.5)};
			bkg.set_walls({wall});
		}
		return bkg;
	}

	// The SOLEDGE3X test mesh's geometry, wall mask and wall outline in a
	// background, without any plasma (which surfaces don't need)
	Background::Background s3x_background(const S3XTest::TestMesh& m,
		const std::string& path)
	{
		S3XTest::write_mesh(path, m);
		const auto mesh {S3X::read_mesh(path, m.nphi)};
		Background::Background bkg {};
		std::vector<BkgFPType> times {0.0, 1.0};
		bkg.move_into_times(times);
		auto zones {mesh.zones};
		bkg.move_into_zones(zones);
		S3X::calc_geometry(mesh, bkg);
		S3X::load_chi(mesh, bkg);
		bkg.set_walls(mesh.walls);
		std::filesystem::remove(path);
		return bkg;
	}

	std::string s3x_mesh_path()
	{
		return (std::filesystem::temp_directory_path()
			/ ("flan_test_sputtering_mesh_" + std::to_string(getpid())
			+ ".h5")).string();
	}
}

// Along a 45 degree wall, every staircase face gets the wall's normal and
// 1 / sqrt(2) of its own area, and is right on the outline
TEST(SputteringSurface, DiagonalWall)
{
	const auto bkg {diagonal_background()};
	const Sputtering::Surface surface(bkg);

	// Plasma cells on the staircase are (i, 9 - i). Each has a wall cell
	// above it (ymax) except the one at the top edge, and to its right
	// (xmax) except the one at the right edge: 9 + 9 faces. The lone wall
	// cell adds 4 more. Every one is in each of the 4 toroidal cells.
	ASSERT_EQ(surface.size(), static_cast<std::size_t>((18 + 4) * diag_nz));
	EXPECT_EQ(surface.num_fallback(), 4 * diag_nz);

	const double s {1.0 / std::sqrt(2.0)};
	double total {};
	for (const auto& e : surface.elements())
	{
		total += e.area;
		if (e.fallback) continue;

		// On the staircase, a face of plasma cell (i, 9 - i)
		EXPECT_EQ(e.i + e.j, diag_n - 1);
		EXPECT_TRUE(e.face == Zones::xmax || e.face == Zones::ymax);

		// Its own area is 1 cm x R at its middle x 0.1 rad
		const double R_mid {e.face == Zones::xmax ? diag_R(e.i + 1)
			: diag_R(e.i + 0.5)};
		EXPECT_NEAR(e.R_mid, R_mid, 1e-12);
		EXPECT_NEAR(e.area_face, 0.01 * R_mid * diag_dphi, 1e-14);

		// The wall's normal, into the plasma (down and to the left), and
		// the area projected onto it
		EXPECT_NEAR(e.n_R, -s, 1e-12);
		EXPECT_NEAR(e.n_Z, -s, 1e-12);
		EXPECT_NEAR(e.area, s * e.area_face, 1e-14);
		EXPECT_NEAR(e.wall_dist, 0.0, 1e-12);
		EXPECT_DOUBLE_EQ(e.phi_max - e.phi_min, diag_dphi);
	}
	EXPECT_NEAR(surface.total_area(), total, 1e-12);
}

// The faces around a wall cell far from the outline fall back to their own
// normal and area. So does every face when there's no outline at all.
TEST(SputteringSurface, Fallback)
{
	const Sputtering::Surface surface(diagonal_background());
	int found {};
	for (const auto& e : surface.elements())
	{
		if (!e.fallback) continue;
		++found;

		// One of the four plasma cells around wall cell (2, 2), with the
		// face toward it, and a normal pointing back into the plasma
		const bool left {e.i == 1 && e.j == 2 && e.face == Zones::xmax};
		const bool right {e.i == 3 && e.j == 2 && e.face == Zones::xmin};
		const bool below {e.i == 2 && e.j == 1 && e.face == Zones::ymax};
		const bool above {e.i == 2 && e.j == 3 && e.face == Zones::ymin};
		EXPECT_TRUE(left || right || below || above);
		EXPECT_NEAR(e.n_R, left ? -1.0 : right ? 1.0 : 0.0, 1e-12);
		EXPECT_NEAR(e.n_Z, below ? -1.0 : above ? 1.0 : 0.0, 1e-12);
		EXPECT_DOUBLE_EQ(e.area, e.area_face);

		// The middle of the face is (x, y) in cell units, and the outline
		// x + y = 10.5 is (10.5 - x - y) / sqrt(2) cells away. That's more
		// than the 2 cell widths (2 cm) allowed.
		const double x {(e.R_mid - 1.0) / 0.01};
		const double y {e.Z_mid / 0.01};
		EXPECT_NEAR(e.wall_dist, 0.01 * (10.5 - x - y) / std::sqrt(2.0),
			1e-12);
		EXPECT_GT(e.wall_dist, 0.02);
	}
	EXPECT_EQ(found, 4 * diag_nz);

	const Sputtering::Surface no_walls(diagonal_background(false));
	EXPECT_EQ(no_walls.size(), surface.size());
	EXPECT_EQ(no_walls.num_fallback(), static_cast<int>(no_walls.size()));
}

// A larger max_wall_dist lets faces farther from the outline use it
TEST(SputteringSurface, MaxWallDist)
{
	const auto bkg {diagonal_background()};
	EXPECT_EQ(Sputtering::Surface(bkg, 100.0).num_fallback(), 0);
}

// Without a wall mask (e.g., a Gkeyll background) there are no surfaces
TEST(SputteringSurface, NeedsChi)
{
	EXPECT_EXIT(Sputtering::Surface {Background::Background {}},
		testing::ExitedWithCode(EXIT_FAILURE), "wall mask");
}

// On the SOLEDGE3X test mesh, the wall is zone 3's last psi row, and the
// outline is a circle along its inner edge. Every face between it and the
// row below is an element, with the circle's inward normal, and they add
// up to the area of that part of the torus.
TEST(SputteringSurface, S3XCircle)
{
	const S3XTest::TestMesh m {};
	const auto bkg {s3x_background(m, s3x_mesh_path())};
	const Sputtering::Surface surface(bkg);

	ASSERT_EQ(surface.size(), static_cast<std::size_t>(m.ntheta * m.nphi));
	EXPECT_EQ(surface.num_fallback(), 0);

	const double r_wall {m.a0 * S3XTest::wall_r(m)};
	for (const auto& e : surface.elements())
	{
		EXPECT_EQ(e.zone, 2);
		EXPECT_EQ(e.i, m.npsi_sol - 2);
		EXPECT_EQ(e.face, Zones::xmax);

		// The normal points in toward the circle's center, (R0, 0)
		const double dR {e.R_mid - m.R0};
		const double r {std::hypot(dR, e.Z_mid)};
		EXPECT_GT(-(e.n_R * dR + e.n_Z * e.Z_mid) / r, 1.0 - 1e-4);
		EXPECT_NEAR(e.area, e.area_face, 1e-4 * e.area_face);
		EXPECT_LT(e.wall_dist, 1e-3 * r_wall);
	}

	// Half a torus of minor radius r_wall (theta from 0 to pi) has area
	// 2 pi^2 R0 r_wall. The faces are chords of the circle, so they're a
	// little short: 1 - (dtheta / 2)^2 / 6, about 3e-4 here.
	const double expected {2.0 * S3XTest::pi * S3XTest::pi * m.R0 * r_wall};
	EXPECT_NEAR(surface.total_area() / expected, 1.0, 1e-3);
}

// A wall cell across a zone face is found from the plasma cell on the
// other side, using the ghost cells. With zone 3's first psi row also
// wall, zone 1's last row has a face on it (and so does zone 3's second
// row, from the other side). Both are far from the outline, so they fall
// back to their own normals.
TEST(SputteringSurface, S3XAcrossZones)
{
	S3XTest::TestMesh m {};
	m.chi_zone3_first_row = true;
	const auto bkg {s3x_background(m, s3x_mesh_path())};
	const Sputtering::Surface surface(bkg);

	const int per_row {m.ntheta * m.nphi};
	ASSERT_EQ(surface.size(), static_cast<std::size_t>(3 * per_row));
	EXPECT_EQ(surface.num_fallback(), 2 * per_row);

	int zone1 {};
	for (const auto& e : surface.elements())
	{
		if (e.zone != 0) continue;
		++zone1;
		EXPECT_EQ(e.i, m.npsi_core - 1);
		EXPECT_EQ(e.face, Zones::xmax);
		EXPECT_TRUE(e.fallback);

		// Pointing back into zone 1, toward the circle's center
		const double dR {e.R_mid - m.R0};
		EXPECT_LT(e.n_R * dR + e.n_Z * e.Z_mid, 0.0);
	}
	EXPECT_EQ(zone1, per_row);
}
