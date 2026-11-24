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
* - The sputtered source on those surfaces (Sputtering::Source): the sheath
*   energy, the effective yield and the ion flux, on the 45 degree wall
*   with a uniform plasma, where they're known exactly.
* - Starting sputtered atoms from the source: the Thompson energies, the
*   cosine directions, the energy of the sputtering ion, and where, when
*   and how fast sample_birth starts them.
*/

#include <cmath>
#include <array>
#include <omp.h>
#include <cstdlib>
#include <sstream>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "background.h"
#include "boundary.h"
#include "constants.h"
#include "flan_types.h"
#include "hdf5.h"
#include "hdf5_utils.h"
#include "options.h"
#include "pcg32.h"
#include "read_bkg.h"
#include "read_s3xe.h"
#include "s3x_test_files.h"
#include "sheath.h"
#include "slots.h"
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
		double projectile_Z {1.0};
		double projectile_mass_amu {2.014};

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

		// A second table, for self-sputtering
		std::string m_self_path {(std::filesystem::temp_directory_path()
			/ ("flan_test_self_sputtering_" + std::to_string(getpid())
			+ ".h5")).string()};

		void write(const Table& t) {write(t, m_path);}

		void write(const Table& t, const std::string& path)
		{
			const hid_t file {H5Fcreate(path.c_str(), H5F_ACC_TRUNC,
				H5P_DEFAULT, H5P_DEFAULT)};
			write_doubles(file, "/energy", {t.energy.size()}, t.energy);
			write_doubles(file, "/angle", {t.angle.size()}, t.angle);
			write_doubles(file, "/yield", t.yield_shape.empty()
				? std::vector<hsize_t> {t.energy.size(), t.angle.size()}
				: t.yield_shape, t.yield);
			write_doubles(file, "/projectile/Z", {}, {t.projectile_Z});
			write_doubles(file, "/projectile/mass_amu", {},
				{t.projectile_mass_amu});
			write_doubles(file, "/target/Z", {}, {74.0});
			write_doubles(file, "/target/mass_amu", {}, {183.84});
			if (t.target_Es) write_doubles(file, "/target/Es", {}, {8.79});
			H5Fclose(file);
		}

		void TearDown() override
		{
			std::filesystem::remove(m_path);
			std::filesystem::remove(m_self_path);
		}
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

	// With chi = false there's no wall mask, and with wall_faces = true
	// the zone's ymin and xmax faces are wall faces (as targets are on a
	// grid fitted to them) and its z faces periodic, with z the toroidal
	// angle.
	Background::Background diagonal_background(bool walls = true,
		bool chi_mask = true, bool wall_faces = false)
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
		if (wall_faces)
		{
			zones.set_face(0, Zones::ymin, Zones::face_wall);
			zones.set_face(0, Zones::xmax, Zones::face_wall);
			zones.set_face(0, Zones::zmin, 0);
			zones.set_face(0, Zones::zmax, 0);
			zones.set_toroidal(0);
		}
		Background::Background bkg {};
		std::vector<BkgFPType> times {0.0, 1.0};
		bkg.move_into_times(times);
		bkg.move_into_zones(zones);

		// Cell centers and corners, in Cartesian coordinates, with
		// (R, z, Z) right handed (the background's default phi_sign = +1)
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
		if (chi_mask) bkg.move_into_chi(chi);

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
		Axisym::calc_geometry(S3X::to_axisym(mesh), bkg);
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

		// On the staircase, a face of plasma cell (i, 9 - i), with the wall
		// cell across it to the right or above
		EXPECT_EQ(e.i + e.j, diag_n - 1);
		EXPECT_TRUE(e.face == Zones::xmax || e.face == Zones::ymax);
		EXPECT_EQ(e.wall_zone, 0);
		EXPECT_EQ(e.wall_i, e.i + (e.face == Zones::xmax ? 1 : 0));
		EXPECT_EQ(e.wall_j, e.j + (e.face == Zones::ymax ? 1 : 0));

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

// Without a wall mask or wall faces on a toroidal grid (e.g., a Gkeyll
// background) there are no surfaces
TEST(SputteringSurface, NeedsChi)
{
	EXPECT_EXIT(Sputtering::Surface {Background::Background {}},
		testing::ExitedWithCode(EXIT_FAILURE), "wall mask");
}

// On a grid fitted to the targets (no wall mask), every cell face on a
// wall face of a zone is an element, with the face's own normal (into the
// plasma) and area. Here those are the bottom row (ymin, normal +Z) and
// the right column (xmax, normal -R) of the 45 degree wall's grid.
TEST(SputteringSurface, TargetFaces)
{
	const auto bkg {diagonal_background(true, false, true)};
	const Sputtering::Surface surface(bkg);
	EXPECT_TRUE(surface.from_wall_faces());
	ASSERT_EQ(surface.size(), static_cast<std::size_t>(2 * diag_n
		* diag_nz));
	EXPECT_EQ(surface.num_fallback(), 0);

	double expected_total {};
	for (const auto& e : surface.elements())
	{
		EXPECT_TRUE(e.on_wall_face());
		EXPECT_EQ(e.zone, 0);
		EXPECT_FALSE(e.fallback);
		EXPECT_DOUBLE_EQ(e.area, e.area_face);
		EXPECT_DOUBLE_EQ(e.phi_min, e.k * diag_dphi);
		EXPECT_DOUBLE_EQ(e.phi_max - e.phi_min, diag_dphi);
		ASSERT_TRUE(e.face == Zones::ymin || e.face == Zones::xmax);
		if (e.face == Zones::ymin)
		{
			EXPECT_EQ(e.j, 0);
			EXPECT_NEAR(e.n_R, 0.0, 1e-12);
			EXPECT_NEAR(e.n_Z, 1.0, 1e-12);
			EXPECT_NEAR(e.Z_mid, 0.0, 1e-14);
			EXPECT_NEAR(e.R_mid, diag_R(e.i + 0.5), 1e-12);
		}
		else
		{
			EXPECT_EQ(e.i, diag_n - 1);
			EXPECT_NEAR(e.n_R, -1.0, 1e-12);
			EXPECT_NEAR(e.n_Z, 0.0, 1e-12);
			EXPECT_NEAR(e.R_mid, diag_R(diag_n), 1e-12);
			EXPECT_NEAR(e.Z_mid, diag_Z(e.j + 0.5), 1e-14);
		}
		const double area {0.01 * e.R_mid * diag_dphi};
		EXPECT_NEAR(e.area, area, 1e-14);
		expected_total += area;
	}
	EXPECT_NEAR(surface.total_area(), expected_total, 1e-12);

	// Only toroidal zones: without z being phi there's no surface
	auto flat {diagonal_background(true, false, true)};
	Zones::Zones zones(2);
	const auto& z0 {flat.get_zones()};
	zones.add_zone(z0.get_x(0), z0.get_y(0), z0.get_z(0), z0.get_grid_x(0),
		z0.get_grid_y(0), z0.get_grid_z(0));
	zones.set_face(0, Zones::ymin, Zones::face_wall);
	flat.move_into_zones(zones);
	EXPECT_EXIT(Sputtering::Surface {flat},
		testing::ExitedWithCode(EXIT_FAILURE), "wall mask");
}

// With a wall mask the wall faces are inside the material, so they're not
// elements: the staircase is the same with or without them
TEST(SputteringSurface, StaircaseIgnoresWallFaces)
{
	const Sputtering::Surface plain(diagonal_background());
	const Sputtering::Surface with(diagonal_background(true, true, true));
	EXPECT_FALSE(with.from_wall_faces());
	ASSERT_EQ(with.size(), plain.size());
	for (std::size_t e {}; e < plain.size(); ++e)
	{
		EXPECT_FALSE(with.elements()[e].on_wall_face());
		EXPECT_EQ(with.elements()[e].i, plain.elements()[e].i);
		EXPECT_EQ(with.elements()[e].j, plain.elements()[e].j);
		EXPECT_EQ(with.elements()[e].face, plain.elements()[e].face);
	}
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

		// The wall cell is across the zone face, in zone 3's first row
		EXPECT_EQ(e.wall_zone, 2);
		EXPECT_EQ(e.wall_i, 0);
		EXPECT_EQ(e.wall_j, e.j);

		// Pointing back into zone 1, toward the circle's center
		const double dR {e.R_mid - m.R0};
		EXPECT_LT(e.n_R * dR + e.n_Z * e.Z_mid, 0.0);
	}
	EXPECT_EQ(zone1, per_row);
}


namespace
{
	// Electron and ion (D) masses for the source tests (amu)
	constexpr double me_amu {5.4858e-4};
	constexpr double mi_amu {2.014};

	// The uniform plasma put on the 45 degree wall background for the
	// source tests, the same in both time frames, in Cartesian components
	struct Plasma
	{
		double n {1e19};		// Density (m^-3)
		double te {20.0};		// Electron temperature (eV)
		double ti {20.0};		// Ion temperature (eV)
		std::array<double, 3> B {0.0, 0.0, 2.0};	// Field (T)
		std::array<double, 3> u {0.0, 0.0, 1e4};	// Flow (m/s)
		std::array<double, 3> E {0.0, 0.0, 0.0};	// Electric field (V/m)
	};

	// The 45 degree wall background with a uniform plasma. B is along Z by
	// default and the flow is up, toward the wall (which is up and to the
	// right of the plasma), with no electric field.
	Background::Background diagonal_with_plasma(const Plasma& p,
		Background::Background bkg = diagonal_background())
	{
		bkg.set_elec_mass_amu(me_amu);
		bkg.set_ion_mass_amu(mi_amu);
		auto field = [&](double value)
		{
			Vectors::ZonedVector<BkgFPType> v(2, {diag_n}, {diag_n},
				{diag_nz}, 2);
			for (auto& x : v.get_data()) x = value;
			return v;
		};
		auto ne {field(p.n)}, te {field(p.te)}, ti {field(p.ti)};
		auto bX {field(p.B[0])}, bY {field(p.B[1])}, bZ {field(p.B[2])};
		auto uX {field(p.u[0])}, uY {field(p.u[1])}, uZ {field(p.u[2])};
		auto eX {field(p.E[0])}, eY {field(p.E[1])}, eZ {field(p.E[2])};
		bkg.move_into_ne(ne);
		bkg.move_into_te(te);
		bkg.move_into_ti(ti);
		bkg.move_into_bX(bX);
		bkg.move_into_bY(bY);
		bkg.move_into_bZ(bZ);
		bkg.move_into_uX(uX);
		bkg.move_into_uY(uY);
		bkg.move_into_uZ(uZ);
		bkg.move_into_eX(eX);
		bkg.move_into_eY(eY);
		bkg.move_into_eZ(eZ);
		return bkg;
	}

	// The element's normal in Cartesian coordinates at the middle of its
	// phi range, the same way Source does it, for a background with
	// toroidal handedness phi_sign (Background::get_phi_sign)
	std::array<double, 3> cartesian_normal(
		const Sputtering::SurfaceElement& e, const double phi_sign)
	{
		const double phi {0.5 * (e.phi_min + e.phi_max)};
		return {e.n_R * std::cos(phi), e.n_R * phi_sign * std::sin(phi),
			e.n_Z};
	}

	// A table with yield c at every energy and angle
	Table constant_table(double c)
	{
		Table t {};
		t.energy = {1.0, 10.0, 100.0, 1e5};
		t.yield.assign(t.energy.size() * t.angle.size(), c);
		return t;
	}

	// A table with yield 0 up to e0 and c from just above it (a step)
	Table step_table(double e0, double c)
	{
		Table t {};
		t.energy = {1.0, e0, e0 * (1.0 + 1e-12), 1e5};
		t.yield.clear();
		for (const double e : t.energy)
		{
			for (std::size_t a {}; a < t.angle.size(); ++a)
				t.yield.push_back(e > e0 ? c : 0.0);
		}
		return t;
	}

	Options::Options source_options()
	{
		Options::Options opts {};
		opts.set_imp_atom_num(74);
		return opts;
	}
}

// The sheath and presheath energy, by hand: for D at Ti = Te it's
// about 3.34 Te
TEST(SputteringSource, SheathEnergy)
{
	const double drop {-0.5 * std::log(2.0 * S3XTest::pi * me_amu / mi_amu
		* 2.0) + 0.5};
	EXPECT_NEAR(drop, 3.3386, 1e-4);
	EXPECT_NEAR(Sputtering::sheath_energy(20.0, 20.0, me_amu, mi_amu, 1),
		20.0 * drop, 1e-12);
	EXPECT_NEAR(Sputtering::sheath_energy(20.0, 20.0, me_amu, mi_amu, 2),
		40.0 * drop, 1e-12);

	// With cold ions, the drop is a little bigger
	const double cold {-0.5 * std::log(2.0 * S3XTest::pi * me_amu / mi_amu)
		+ 0.5};
	EXPECT_NEAR(Sputtering::sheath_energy(10.0, 0.0, me_amu, mi_amu, 1),
		10.0 * cold, 1e-12);
}

// A constant yield averages to itself. A step from 0 to c at e0 averages
// to c times the fraction of a Gamma(2, Ti) distribution above
// x0 = (e0 - E_sheath) / Ti, which is (1 + x0) e^-x0. The step is put on a
// boundary between the quadrature's points, x0 = 1 or 2, so the midpoint
// rule is good to about dx^2.
TEST_F(SputteringTable, EffectiveYield)
{
	write(constant_table(0.02));
	const Sputtering::YieldTable constant(m_path);
	EXPECT_NEAR(Sputtering::effective_yield(constant, 60.0, 20.0, 45.0)
		.yield, 0.02, 1e-14);

	for (const double x0 : {1.0, 2.0})
	{
		const double e_sheath {60.0};
		const double ti {20.0};
		write(step_table(e_sheath + x0 * ti, 0.05));
		const Sputtering::YieldTable step(m_path);
		const auto y {Sputtering::effective_yield(step, e_sheath, ti,
			30.0)};
		EXPECT_NEAR(y.yield, 0.05 * (1.0 + x0) * std::exp(-x0), 2e-3 * 0.05)
			<< "x0 = " << x0;
		EXPECT_EQ(y.frac_above, 0.0);
	}

	// Everything below the threshold: no yield
	write(step_table(1e4, 0.05));
	EXPECT_NEAR(Sputtering::effective_yield(Sputtering::YieldTable {m_path},
		60.0, 20.0, 30.0).yield, 0.0, 1e-12);
}

// Parallel flux onto the 45 degree wall, with B along Z and the flow up:
// n u |b . n| = n u / sqrt(2), and the rest of the source follows from it
TEST_F(SputteringTable, SourceParallel)
{
	write(constant_table(0.02));
	const Sputtering::YieldTable table(m_path);
	const Plasma p {};
	const auto bkg {diagonal_with_plasma(p)};
	const Sputtering::Surface surface(bkg);
	const Sputtering::Source source(surface, bkg, table, source_options());

	ASSERT_EQ(source.size(), surface.size());
	ASSERT_EQ(source.nframes(), 2);
	const double e_sheath {Sputtering::sheath_energy(p.te, p.ti, me_amu,
		mi_amu, 1)};
	for (int t {}; t < 2; ++t)
	{
		double total {};
		for (std::size_t e {}; e < source.size(); ++e)
		{
			const auto& el {surface.elements()[e]};
			const std::size_t idx {t * source.size() + e};
			const auto n {cartesian_normal(el, bkg.get_phi_sign())};

			// u_par = u_z, b_hat = Z_hat, so the parallel velocity into
			// the wall is -u_z n_z
			const double expected {p.n * std::max(0.0, -p.u[2] * n[2])};
			EXPECT_NEAR(source.flux_par()[idx], -p.n * p.u[2] * n[2],
				1e-6 * p.n * p.u[2]);
			EXPECT_NEAR(source.flux()[idx], expected, 1e-6 * p.n * p.u[2]);
			EXPECT_NEAR(source.flux_exb()[idx], 0.0, 1e-6);
			EXPECT_NEAR(source.energy()[idx], e_sheath + 2.0 * p.ti, 1e-10);
			EXPECT_NEAR(source.yield()[idx], 0.02, 1e-14);
			EXPECT_NEAR(source.rate()[idx], expected * 0.02 * el.area,
				1e-6 * p.n * p.u[2] * 0.02 * el.area);
			total += source.rate()[idx];

			// On the outline, n_Z = -1 / sqrt(2), so it's n u / sqrt(2)
			if (!el.fallback)
			{
				EXPECT_NEAR(source.flux()[idx], p.n * p.u[2]
					/ std::sqrt(2.0), 1e-6 * p.n * p.u[2]);
			}
		}
		EXPECT_NEAR(source.total_rate()[t], total, 1e-9 * total);
	}
	EXPECT_NEAR(source.mean_total_rate(), source.total_rate()[0],
		1e-9 * source.total_rate()[0]);
	EXPECT_EQ(source.num_above_table(), 0);
}

// The E x B drift adds to the flux onto the wall. With B along Z and E
// along Y, v_E = E / B along X. "parallel" leaves it out. Flow away from
// the wall gives no flux.
TEST_F(SputteringTable, SourceExB)
{
	write(constant_table(0.02));
	const Sputtering::YieldTable table(m_path);
	Plasma p {};
	p.E = {0.0, 5e3, 0.0};
	const double v_E {p.E[1] / p.B[2]};
	const auto bkg {diagonal_with_plasma(p)};
	const Sputtering::Surface surface(bkg);
	const Sputtering::Source source(surface, bkg, table, source_options());

	auto opts_par {source_options()};
	opts_par.set_sput_flux_opt("parallel");
	const Sputtering::Source par(surface, bkg, table, opts_par);

	for (std::size_t e {}; e < source.size(); ++e)
	{
		const auto n {cartesian_normal(surface.elements()[e],
			bkg.get_phi_sign())};
		const double f_par {-p.n * p.u[2] * n[2]};
		const double f_exb {-p.n * v_E * n[0]};
		const double scale {p.n * p.u[2]};
		EXPECT_NEAR(source.flux_exb()[e], f_exb, 1e-6 * scale);
		EXPECT_NEAR(source.flux()[e], std::max(0.0, f_par + f_exb),
			1e-6 * scale);
		EXPECT_NEAR(par.flux()[e], std::max(0.0, f_par), 1e-6 * scale);
		EXPECT_NEAR(par.flux_exb()[e], f_exb, 1e-6 * scale);
	}

	// Flow down, away from the wall: no parallel flux into it
	Plasma away {};
	away.u = {0.0, 0.0, -1e4};
	const auto bkg_away {diagonal_with_plasma(away)};
	const Sputtering::Source none(Sputtering::Surface {bkg_away}, bkg_away,
		table, source_options());
	for (std::size_t e {}; e < none.size(); ++e)
	{
		const auto n {cartesian_normal(none.surface().elements()[e],
			bkg_away.get_phi_sign())};
		if (n[2] < 0.0) EXPECT_EQ(none.flux()[e], 0.0);
	}
}

// With the Bohm floor, a slow parallel flow onto the wall is replaced by
// the sound speed, sqrt(e (Te + Ti) / m_i)
TEST_F(SputteringTable, SourceBohmFloor)
{
	write(constant_table(0.02));
	const Sputtering::YieldTable table(m_path);
	Plasma p {};
	p.u = {0.0, 0.0, 10.0};
	const auto bkg {diagonal_with_plasma(p)};
	const Sputtering::Surface surface(bkg);
	auto opts {source_options()};
	opts.set_sput_bohm_floor("on");
	const Sputtering::Source source(surface, bkg, table, opts);

	const double c_s {std::sqrt(1.602e-19 * (p.te + p.ti)
		/ (mi_amu * 1.661e-27))};
	for (std::size_t e {}; e < source.size(); ++e)
	{
		const auto n {cartesian_normal(surface.elements()[e],
			bkg.get_phi_sign())};
		EXPECT_NEAR(source.flux_par()[e], p.n * std::max(-p.u[2] * n[2],
			c_s * std::abs(n[2])), 1e-6 * p.n * c_s);
	}
}

// A table for the wrong impurity or projectile, or a bad angle, stops the
// run. A projectile mass far from the background's, or energies past the
// end of the table, are reported.
TEST_F(SputteringTable, SourceChecks)
{
	const auto bkg {diagonal_with_plasma(Plasma {})};
	const Sputtering::Surface surface(bkg);
	write(constant_table(0.02));
	const Sputtering::YieldTable table(m_path);

	auto carbon {source_options()};
	carbon.set_imp_atom_num(6);
	EXPECT_EXIT(Sputtering::Source(surface, bkg, table, carbon),
		testing::ExitedWithCode(EXIT_FAILURE), "imp_atom_num");

	auto angle {source_options()};
	angle.set_sput_angle_val(95.0);
	EXPECT_EXIT(Sputtering::Source(surface, bkg, table, angle),
		testing::ExitedWithCode(EXIT_FAILURE), "sput_angle_val");

	Table helium {constant_table(0.02)};
	helium.projectile_Z = 2.0;
	write(helium);
	EXPECT_EXIT(Sputtering::Source(surface, bkg,
		Sputtering::YieldTable {m_path}, source_options()),
		testing::ExitedWithCode(EXIT_FAILURE), "Z = 2");

	// Heavier background ions than the table's projectile
	auto heavy {diagonal_with_plasma(Plasma {})};
	heavy.set_ion_mass_amu(4.0);
	const Sputtering::Source mass(surface, heavy, table, source_options());
	std::ostringstream out {};
	mass.print_summary(out);
	EXPECT_NE(out.str().find("more than 5% off"), std::string::npos);

	// A table that ends at 50 eV, below every impact energy
	Table low {constant_table(0.02)};
	low.energy = {1.0, 10.0, 20.0, 50.0};
	write(low);
	const Sputtering::Source above(surface, bkg,
		Sputtering::YieldTable {m_path}, source_options());
	EXPECT_EQ(above.num_above_table(),
		static_cast<int>(2 * surface.size()));
}

// Thompson energies come out exactly as the inverse of the cut-off
// cumulative distribution: never above e_max, the median where F(E) is half
// of F(e_max), and the mean that of f(E) ~ E / (E + e_b)^3 on [0, e_max]
TEST(SputteringBirth, Thompson)
{
	const double e_b {8.79};
	const double e_max {20.0};
	auto F = [&](double e) {return e * e / ((e + e_b) * (e + e_b));};

	const int n {200000};
	double mean {};
	for (int k {}; k < n; ++k)
	{
		const double e {Sputtering::sample_thompson(e_b, e_max,
			(k + 0.5) / n)};
		EXPECT_LE(e, e_max);
		EXPECT_GE(e, 0.0);
		mean += e / n;
	}
	EXPECT_NEAR(F(Sputtering::sample_thompson(e_b, e_max, 0.5)),
		0.5 * F(e_max), 1e-12);

	// Mean of f on [0, e_max], numerically
	double num {}, den {};
	const int m {200000};
	for (int k {}; k < m; ++k)
	{
		const double e {(k + 0.5) * e_max / m};
		const double f {e / std::pow(e + e_b, 3)};
		num += e * f;
		den += f;
	}
	EXPECT_NEAR(mean, num / den, 1e-3 * num / den);

	// No room for any energy: the binding energy
	EXPECT_EQ(Sputtering::sample_thompson(e_b, 0.0, 0.3), e_b);
}

// Cosine directions are unit vectors on the normal's side, with cos(theta)
// = sqrt(u1), so cos^2(theta) averages 1/2, and no bias sideways
TEST(SputteringBirth, CosineDirection)
{
	std::array<double, 3> n {1.0, 2.0, -0.5};
	const double len {std::sqrt(1.0 + 4.0 + 0.25)};
	for (double& c : n) c /= len;

	const int m {200};
	double cos2 {};
	std::array<double, 3> sum {};
	for (int a {}; a < m; ++a)
	for (int b {}; b < m; ++b)
	{
		const double u1 {(a + 0.5) / m};
		const auto d {Sputtering::cosine_direction(n, u1, (b + 0.5) / m)};
		const double dn {d[0] * n[0] + d[1] * n[1] + d[2] * n[2]};
		EXPECT_NEAR(d[0] * d[0] + d[1] * d[1] + d[2] * d[2], 1.0, 1e-12);
		EXPECT_NEAR(dn, std::sqrt(u1), 1e-12);
		cos2 += dn * dn / (m * m);
		for (int c {}; c < 3; ++c) sum[c] += (d[c] - dn * n[c]) / (m * m);
	}
	EXPECT_NEAR(cos2, 0.5, 1e-4);
	for (int c {}; c < 3; ++c) EXPECT_NEAR(sum[c], 0.0, 1e-12);
}

// The sputtering ions' energies follow Y(E) p(eps). With a step table (0
// below E0, c above), none is below E0, and x = eps / Ti has the mean of
// x e^-x above x0, (x0^2 + 2 x0 + 2) / (x0 + 1). With no yield anywhere,
// the mean impact energy comes back.
TEST_F(SputteringTable, ImpactEnergy)
{
	const double e_sheath {60.0};
	const double ti {20.0};
	const double x0 {1.0};
	write(step_table(e_sheath + x0 * ti, 0.05));
	const Sputtering::YieldTable table(m_path);

	const int n {20000};
	double mean_x {};
	for (int k {}; k < n; ++k)
	{
		const double e {Sputtering::sample_impact_energy(table, e_sheath, ti,
			30.0, (k + 0.5) / n)};
		EXPECT_GT(e, e_sheath + x0 * ti);
		mean_x += (e - e_sheath) / ti / n;
	}
	EXPECT_NEAR(mean_x, (x0 * x0 + 2.0 * x0 + 2.0) / (x0 + 1.0), 5e-3);

	write(step_table(1e4, 0.05));
	EXPECT_DOUBLE_EQ(Sputtering::sample_impact_energy(
		Sputtering::YieldTable {m_path}, e_sheath, ti, 30.0, 0.7),
		e_sheath + 2.0 * ti);
}

// On the 45 degree wall: births pick frames and elements in proportion to
// the atoms they sputter, and start on the right face, at the right time,
// moving into the plasma with an energy below the cut-off
TEST_F(SputteringTable, SampleBirth)
{
	write(constant_table(0.02));
	const Sputtering::YieldTable table(m_path);
	const auto bkg {diagonal_with_plasma(Plasma {})};
	const Sputtering::Surface surface(bkg);
	const Sputtering::Source source(surface, bkg, table, source_options());

	// Two frames at t = 0 and 1, each standing for half the time
	EXPECT_EQ(source.frame_lo(), (std::vector<double> {0.0, 0.5}));
	EXPECT_EQ(source.frame_hi(), (std::vector<double> {0.5, 1.0}));
	EXPECT_NEAR(source.source_rate(), source.mean_total_rate(),
		1e-12 * source.source_rate());

	const std::size_t n_elem {source.size()};
	const double total {0.5 * (source.total_rate()[0]
		+ source.total_rate()[1])};
	const double lambda {4.0 * 2.014 * 183.84 / std::pow(2.014 + 183.84, 2)};
	const double m_w {183.84 * 1.661e-27};

	pcg32 rng(42, 7);
	const int n {50000};
	std::vector<int> counts(2 * n_elem, 0);
	for (int k {}; k < n; ++k)
	{
		const auto b {source.sample_birth(rng)};
		++counts[b.frame * n_elem + b.element];
		const auto& el {surface.elements()[b.element]};

		// Time in the frame's window, phi in the element's range
		EXPECT_GE(b.t, source.frame_lo()[b.frame]);
		EXPECT_LE(b.t, source.frame_hi()[b.frame]);
		EXPECT_GE(b.z, el.phi_min);
		EXPECT_LT(b.z, el.phi_max);

		// In the plasma cell, a thousandth of a cell from the face (x and
		// y are cell indices here)
		EXPECT_EQ(b.zone, 0);
		EXPECT_GE(b.x, el.i);
		EXPECT_LE(b.x, el.i + 1);
		EXPECT_GE(b.y, el.j);
		EXPECT_LE(b.y, el.j + 1);
		if (el.face == Zones::xmax) EXPECT_NEAR(b.x, el.i + 0.999, 1e-12);
		if (el.face == Zones::ymax) EXPECT_NEAR(b.y, el.j + 0.999, 1e-12);
		if (el.face == Zones::xmin) EXPECT_NEAR(b.x, el.i + 0.001, 1e-12);
		if (el.face == Zones::ymin) EXPECT_NEAR(b.y, el.j + 0.001, 1e-12);

		// Speed from the atom's energy, below the cut-off, moving into the
		// plasma side of the element's normal at its phi
		EXPECT_LE(b.e_atom, lambda * b.e_impact * (1.0 + 1e-12));
		const double v2 {b.vX * b.vX + b.vY * b.vY + b.vZ * b.vZ};
		EXPECT_NEAR(0.5 * m_w * v2 / 1.602e-19, b.e_atom,
			1e-9 * b.e_atom + 1e-12);
		const double n_x {el.n_R * std::cos(b.z)};
		const double n_y {el.n_R * bkg.get_phi_sign() * std::sin(b.z)};
		EXPECT_GE(b.vX * n_x + b.vY * n_y + b.vZ * el.n_Z, 0.0);
	}

	// Chi-square of the counts against rate x window / total
	double chi2 {};
	int bins {};
	for (int t {}; t < 2; ++t)
	for (std::size_t e {}; e < n_elem; ++e)
	{
		const double expected {n * source.rate()[t * n_elem + e] * 0.5
			/ total};
		if (expected <= 0.0)
		{
			EXPECT_EQ(counts[t * n_elem + e], 0);
			continue;
		}
		const double d {counts[t * n_elem + e] - expected};
		chi2 += d * d / expected;
		++bins;
	}
	EXPECT_LT(chi2 / (bins - 1), 1.5) << "chi2 = " << chi2 << " for "
		<< bins << " bins";
}

// On target faces, with B along Z and the flow down onto the bottom
// target: the bottom row gets n u (b . n = 1) and the right column, which
// the field runs along, nothing. Births start just above the bottom
// faces, moving up into the plasma.
TEST_F(SputteringTable, SourceOnTargets)
{
	write(constant_table(0.02));
	const Sputtering::YieldTable table(m_path);
	Plasma p {};
	p.u = {0.0, 0.0, -1e4};
	const auto bkg {diagonal_with_plasma(p, diagonal_background(true, false,
		true))};
	const Sputtering::Surface surface(bkg);
	const Sputtering::Source source(surface, bkg, table, source_options());

	double total {};
	for (std::size_t e {}; e < source.size(); ++e)
	{
		const auto& el {surface.elements()[e]};
		const double flux {el.face == Zones::ymin ? p.n * 1e4 : 0.0};
		EXPECT_NEAR(source.flux()[e], flux, 1e-6 * p.n * 1e4);
		EXPECT_NEAR(source.rate()[e], flux * 0.02 * el.area,
			1e-6 * p.n * 1e4 * 0.02 * el.area);
		if (el.face == Zones::ymin) total += el.area;
	}
	EXPECT_NEAR(source.total_rate()[0], p.n * 1e4 * 0.02 * total,
		1e-9 * source.total_rate()[0]);

	pcg32 rng(5, 3);
	for (int n {}; n < 2000; ++n)
	{
		const auto b {source.sample_birth(rng)};
		const auto& el {surface.elements()[b.element]};
		ASSERT_EQ(el.face, Zones::ymin);
		EXPECT_NEAR(b.y, 0.001, 1e-12);
		EXPECT_GE(b.x, el.i);
		EXPECT_LE(b.x, el.i + 1);
		EXPECT_GE(b.vZ, 0.0);
	}
}

// Filling the slots from a source gives neutral particles where the source
// put them, each standing for source_rate / imp_num atoms/s
TEST_F(SputteringTable, FillSlots)
{
	write(constant_table(0.02));
	const Sputtering::YieldTable table(m_path);
	const auto bkg {diagonal_with_plasma(Plasma {})};
	auto opts {source_options()};
	opts.set_imp_num(50);
	opts.set_imp_init_charge(3);
	const Sputtering::Source source(Sputtering::Surface {bkg}, bkg, table,
		opts);

	Slots::Slots slots(50);
	for (int i {}; i < slots.N(); ++i) slots.set_state(i, 1);
	std::vector<pcg32> rngs {};
	for (int i {}; i < omp_get_max_threads(); ++i) rngs.emplace_back(3, i);
	int rem {50};
	int alive {};
	Slots::fill_slots_cpu(slots, rem, alive, bkg, opts, rngs, &source);

	EXPECT_EQ(alive, 50);
	EXPECT_EQ(rem, 0);
	for (int i {}; i < slots.N(); ++i)
	{
		EXPECT_EQ(slots.q()[i], 0);
		EXPECT_DOUBLE_EQ(slots.weight()[i], source.source_rate() / 50);
		EXPECT_EQ(slots.zone()[i], 0);
		EXPECT_EQ(slots.xidx()[i], static_cast<int>(slots.x()[i]));
		EXPECT_EQ(slots.yidx()[i], static_cast<int>(slots.y()[i]));
	}
}

// With nothing sputtered, there's nothing to start
TEST_F(SputteringTable, SampleBirthNoSource)
{
	write(constant_table(0.0));
	const auto bkg {diagonal_with_plasma(Plasma {})};
	const Sputtering::Source source(Sputtering::Surface {bkg}, bkg,
		Sputtering::YieldTable {m_path}, source_options());
	pcg32 rng(1, 1);
	EXPECT_EXIT(source.sample_birth(rng),
		testing::ExitedWithCode(EXIT_FAILURE), "nothing is sputtered");
}

// An atom absorbed in a wall cell goes on the element between that wall
// cell and the cell it came from. A wall cell on the 45 degree staircase
// borders two plasma cells, so which element depends on where it came
// from. A wall cell with no element is unassigned. The time since birth
// and charge are binned, and wall faces (with no elements on them here)
// counted on their own.
TEST(SputteringDeposition, Tallies)
{
	const auto bkg {diagonal_background()};
	const Sputtering::Surface surface(bkg);
	Sputtering::Deposition dep(surface, 74);

	// Wall cell (5, 5) borders plasma cells (4, 5) (its xmin side) and
	// (5, 4) (its ymin side), in toroidal cell 1
	auto find = [&](int i, int j, Zones::Face face)
	{
		for (std::size_t e {}; e < surface.size(); ++e)
		{
			const auto& el {surface.elements()[e]};
			if (el.i == i && el.j == j && el.k == 1 && el.face == face)
				return e;
		}
		return surface.size();
	};
	const std::size_t from_left {find(4, 5, Zones::xmax)};
	const std::size_t from_below {find(5, 4, Zones::ymax)};
	ASSERT_LT(from_left, surface.size());
	ASSERT_LT(from_below, surface.size());

	dep.add_into_wall_cell(0, 5, 5, 1, 0, 4, 5, 2e-6, 3, 10.0);
	dep.add_into_wall_cell(0, 5, 5, 1, 0, 5, 4, 2e-6, 3, 20.0);
	dep.add_into_wall_cell(0, 9, 9, 1, 0, 8, 9, 5e-4, 0, 40.0);
	dep.add_at_wall_face(bkg.get_zones(), 0, Zones::ymin,
		{1.5, -0.01, 0.15}, 3e-12, 74, 80.0);

	EXPECT_DOUBLE_EQ(dep.redeposited()[from_left], 10.0);
	EXPECT_DOUBLE_EQ(dep.redeposited()[from_below], 20.0);
	EXPECT_DOUBLE_EQ(dep.unassigned(), 40.0);
	EXPECT_DOUBLE_EQ(dep.wall_face(), 80.0);
	EXPECT_DOUBLE_EQ(dep.total(), 150.0);

	// 2e-6 s is in the bin starting at 10^-5.7, 5e-4 s at 10^-3.4, and
	// shorter than 1e-10 s in the first bin
	const auto& edges {dep.lifetime_edges()};
	auto bin_of = [&](double t)
	{
		for (std::size_t b {}; b + 1 < edges.size(); ++b)
			if (edges[b] <= t && t < edges[b + 1]) return b;
		return edges.size();
	};
	EXPECT_DOUBLE_EQ(dep.lifetime()[bin_of(2e-6)], 30.0);
	EXPECT_DOUBLE_EQ(dep.lifetime()[bin_of(5e-4)], 40.0);
	EXPECT_DOUBLE_EQ(dep.lifetime()[0], 80.0);
	EXPECT_DOUBLE_EQ(dep.charge()[3], 30.0);
	EXPECT_DOUBLE_EQ(dep.charge()[0], 40.0);
	EXPECT_DOUBLE_EQ(dep.charge()[74], 80.0);

	std::ostringstream out {};
	dep.print_summary(out, 300.0);
	EXPECT_NE(out.str().find("50% of the source"), std::string::npos);
}

// On target faces, an atom absorbed at a wall face goes on the element at
// the cell it crossed, found from its position along the face and in z
// (wrapped if it's just past the end of the periodic z range). A wall face
// with no elements, or no wall face at all, isn't on an element.
TEST(SputteringDeposition, WallFaces)
{
	const auto bkg {diagonal_background(true, false, true)};
	const Zones::Zones& zones {bkg.get_zones()};
	const Sputtering::Surface surface(bkg);
	Sputtering::Deposition dep(surface, 74);

	auto find = [&](int i, int j, int k, Zones::Face face)
	{
		for (std::size_t e {}; e < surface.size(); ++e)
		{
			const auto& el {surface.elements()[e]};
			if (el.i == i && el.j == j && el.k == k && el.face == face)
				return e;
		}
		return surface.size();
	};
	const std::size_t bottom {find(3, 0, 2, Zones::ymin)};
	const std::size_t wrapped {find(7, 0, 0, Zones::ymin)};
	const std::size_t right {find(diag_n - 1, 6, 1, Zones::xmax)};
	ASSERT_LT(bottom, surface.size());
	ASSERT_LT(wrapped, surface.size());
	ASSERT_LT(right, surface.size());

	dep.add_at_wall_face(zones, 0, Zones::ymin, {3.4, -0.01, 0.25}, 1e-6, 1,
		1.0);
	dep.add_at_wall_face(zones, 0, Zones::ymin,
		{7.2, -0.02, diag_nz * diag_dphi + 0.01}, 1e-6, 1, 2.0);
	dep.add_at_wall_face(zones, 0, Zones::xmax, {diag_n + 0.01, 6.5, 0.15},
		1e-6, 1, 4.0);
	dep.add_at_wall_face(zones, 0, Zones::ymax, {3.5, diag_n + 0.01, 0.15},
		1e-6, 1, 8.0);

	EXPECT_DOUBLE_EQ(dep.redeposited()[bottom], 1.0);
	EXPECT_DOUBLE_EQ(dep.redeposited()[wrapped], 2.0);
	EXPECT_DOUBLE_EQ(dep.redeposited()[right], 4.0);
	EXPECT_DOUBLE_EQ(dep.wall_face(), 8.0);
	EXPECT_DOUBLE_EQ(dep.total(), 15.0);
}

namespace
{
	// Slab options with a wall on one face. The faces default to
	// absorbing, so the one across from the wall isn't periodic.
	Options::Options slab_wall_options(const std::string& wall,
		const double b_angle = 0.0)
	{
		Options::Options opts {source_options()};
		opts.set_bkg_source("slab");
		opts.set_slab_ne(1e19);
		opts.set_slab_te(20.0);
		opts.set_slab_ti(20.0);
		opts.set_slab_bmag(2.0);
		opts.set_slab_tmax(1e-5);
		opts.set_slab_wall(wall);
		opts.set_slab_b_angle(b_angle);
		return opts;
	}
}

// A slab wall at y min is the 8 x 8 cell faces of that side of the 5 cm
// box: normal +Y, each 6.25 mm square, 25 cm^2 in all. It is a Cartesian
// surface, so the normal doesn't turn with z.
TEST(SlabWall, SurfaceYmin)
{
	const auto bkg {Background::read_bkg(slab_wall_options("ymin"))};
	ASSERT_EQ(bkg.get_zones().get_face(0, Zones::ymin), Zones::face_wall);
	EXPECT_TRUE(bkg.get_zones().is_cartesian(0));

	const Sputtering::Surface surface(bkg);
	EXPECT_TRUE(surface.from_wall_faces());
	ASSERT_EQ(surface.size(), 64u);
	const double side {0.05 / 8.0};
	for (const auto& e : surface.elements())
	{
		EXPECT_TRUE(e.cartesian);
		EXPECT_TRUE(e.on_wall_face());
		EXPECT_EQ(e.face, Zones::ymin);
		EXPECT_EQ(e.j, 0);
		EXPECT_NEAR(e.area, side * side, 1e-15);
		EXPECT_NEAR(e.Z_mid, -0.025, 1e-15);	// Y of the face
		EXPECT_NEAR(e.R_mid, (e.i + 0.5) * side, 1e-15);	// X
		EXPECT_NEAR(e.phi_min, -0.025 + e.k * side, 1e-15);	// Z
		for (const double z : {e.phi_min, e.phi_max})
		{
			const auto n {e.normal_xyz(z, 1.0)};
			EXPECT_NEAR(n[0], 0.0, 1e-15);
			EXPECT_NEAR(n[1], 1.0, 1e-15);
			EXPECT_NEAR(n[2], 0.0, 1e-15);
		}
	}
	EXPECT_NEAR(surface.total_area(), 0.05 * 0.05, 1e-14);
}

// At x max the normal is -X, and the surface is the y-z side
TEST(SlabWall, SurfaceXmax)
{
	const auto bkg {Background::read_bkg(slab_wall_options("xmax"))};
	const Sputtering::Surface surface(bkg);
	ASSERT_EQ(surface.size(), 64u);
	for (const auto& e : surface.elements())
	{
		EXPECT_EQ(e.face, Zones::xmax);
		EXPECT_EQ(e.i, 7);
		EXPECT_NEAR(e.R_mid, 0.05, 1e-15);
		const auto n {e.normal_xyz(0.0, 1.0)};
		EXPECT_NEAR(n[0], -1.0, 1e-15);
		EXPECT_NEAR(n[1], 0.0, 1e-15);
	}
	EXPECT_NEAR(surface.total_area(), 0.05 * 0.05, 1e-14);
}

// slab_b_angle turns B from +X toward -Y. With the wall at y min it is the
// angle between B and the wall, and the flow along B runs into the wall.
TEST(SlabWall, FieldAngle)
{
	for (const double deg : {0.0, 3.0, 30.0, 87.0, 90.0})
	{
		auto opts {slab_wall_options("ymin", deg)};
		opts.set_slab_flow(1e4);
		const auto bkg {Background::read_bkg(opts)};
		const double a {deg * S3XTest::pi / 180.0};
		for (const int t : {0, 127})
		{
			const double bX {bkg.get_bX()(0, t, 3, 4, 5)};
			const double bY {bkg.get_bY()(0, t, 3, 4, 5)};
			EXPECT_NEAR(bX, 2.0 * std::cos(a), 1e-6);
			EXPECT_NEAR(bY, -2.0 * std::sin(a), 1e-6);
			EXPECT_NEAR(bkg.get_bZ()(0, t, 3, 4, 5), 0.0, 1e-12);
			EXPECT_NEAR(bkg.get_bmag()(0, t, 3, 4, 5), 2.0, 1e-6);

			// Angle between B and the wall (normal +Y)
			EXPECT_NEAR(std::asin(std::abs(bY) / 2.0), a, 1e-6);
			EXPECT_NEAR(bkg.get_uX()(0, t, 3, 4, 5), 1e4 * std::cos(a), 1e-2);
			EXPECT_NEAR(bkg.get_uY()(0, t, 3, 4, 5), -1e4 * std::sin(a), 1e-2);
		}
	}
}

// The plasma flowing along B at angle a to the wall puts a flux n u sin a
// on it, and atoms are born just off the wall, moving away from it
TEST_F(SputteringTable, SlabWallSource)
{
	write(constant_table(0.02));
	const Sputtering::YieldTable table(m_path);
	const double deg {10.0};
	auto opts {slab_wall_options("ymin", deg)};
	opts.set_slab_flow(1e4);
	const auto bkg {Background::read_bkg(opts)};
	const Sputtering::Surface surface(bkg);
	const Sputtering::Source source(surface, bkg, table, opts);

	const double flux {1e19 * 1e4 * std::sin(deg * S3XTest::pi / 180.0)};
	for (std::size_t e {}; e < source.size(); ++e)
	{
		EXPECT_NEAR(source.flux()[e], flux, 1e-5 * flux);
		EXPECT_NEAR(source.rate()[e], flux * 0.02
			* surface.elements()[e].area, 1e-5 * flux * 0.02
			* surface.elements()[e].area);
	}
	EXPECT_NEAR(source.total_rate()[0], flux * 0.02 * 0.05 * 0.05,
		1e-5 * flux * 0.02 * 0.05 * 0.05);

	pcg32 rng(7, 1);
	const double side {0.05 / 8.0};
	for (int n {}; n < 2000; ++n)
	{
		const auto b {source.sample_birth(rng)};
		const auto& el {surface.elements()[b.element]};
		EXPECT_NEAR(b.y, -0.025 + 1e-3 * side, 1e-12);
		EXPECT_GE(b.x, el.R0 - 1e-15);
		EXPECT_LE(b.x, el.R0 + side + 1e-15);
		EXPECT_GE(b.z, el.phi_min);
		EXPECT_LE(b.z, el.phi_max);
		EXPECT_GE(b.vY, 0.0);
	}
}

// D+ impacts from the sheath model (sput_impact_opt = "sheath"): with a
// constant yield the averaged yield is exactly that, whatever the
// impacts, so the rate is the same as the formula's; the slab's uniform
// plasma is one point of the impact table; the ions hit nearer the
// normal than the field (80 degrees from it here); and each birth
// carries an ion energy and angle.
//
// The mean impact energy is the drop plus what the flux into the sheath
// carries. For a Maxwellian at Ti drifting at c_s along B, weighted by
// the speed along B, that's m c_s^2 / 2 + 2.5 Ti (the speed along B adds
// m c_s^2 / 2 + 1.5 Ti, across it Ti): with c_s^2 = (Te + Ti) e / m,
// 0.5 Te + 3 Ti on top of the floating drop. The formula has 0.5 Te +
// 2 Ti, so the sheath model's is Ti more (only ions moving toward the
// wall are kept, which changes it little at Ti = Te).
TEST_F(SputteringTable, SheathImpacts)
{
	write(constant_table(0.02));
	const Sputtering::YieldTable table(m_path);
	auto opts {slab_wall_options("ymin", 10.0)};
	opts.set_slab_flow(1e4);
	const auto bkg {Background::read_bkg(opts)};
	const Sputtering::Surface surface(bkg);
	const Sputtering::Source formula(surface, bkg, table, opts);
	const Sputtering::Source sheath(surface, bkg, table, opts,
		Sheath::Model::brooks);

	EXPECT_EQ(formula.impacts(), nullptr);
	ASSERT_NE(sheath.impacts(), nullptr);
	EXPECT_EQ(sheath.impacts()->size(), 1u);
	EXPECT_TRUE(sheath.impacts()->built());
	for (std::size_t e {}; e < sheath.size(); ++e)
	{
		EXPECT_NEAR(sheath.yield()[e], 0.02, 1e-12);
		EXPECT_NEAR(sheath.rate()[e], formula.rate()[e],
			1e-12 * formula.rate()[e]);
		EXPECT_LT(sheath.impact_angle()[e], 70.0);
		EXPECT_GT(sheath.impact_angle()[e], 5.0);
		EXPECT_NEAR(formula.impact_angle()[e], 60.0, 1e-12);
		// The formula's is the drop + 0.5 Te + 2 Ti
		const double expected {formula.energy()[e] + 20.0};
		EXPECT_NEAR(sheath.energy()[e], expected, 0.03 * expected);
	}
	EXPECT_NEAR(sheath.source_rate(), formula.source_rate(),
		1e-12 * formula.source_rate());

	// The histograms of the ions hitting the wall add up to the ion flux
	// onto it, and weighted by the yield to the source. Only with the
	// sheath model.
	EXPECT_TRUE(formula.ion_energy().empty());
	auto sum = [](const std::vector<double>& v)
	{
		double total {};
		for (const double x : v) total += x;
		return total;
	};
	const double ions {sheath.source_rate() / 0.02};
	EXPECT_EQ(sheath.ion_energy_edges().size(),
		sheath.ion_energy().size() + 1);
	EXPECT_EQ(sheath.ion_angle().size(), 90u);
	EXPECT_NEAR(sum(sheath.ion_energy()), ions, 1e-9 * ions);
	EXPECT_NEAR(sum(sheath.ion_angle()), ions, 1e-9 * ions);
	EXPECT_NEAR(sum(sheath.sputtering_ion_energy()), sheath.source_rate(),
		1e-9 * sheath.source_rate());
	EXPECT_NEAR(sum(sheath.sputtering_ion_angle()), sheath.source_rate(),
		1e-9 * sheath.source_rate());

	pcg32 rng(5, 1);
	for (int n {}; n < 200; ++n)
	{
		const auto b {sheath.sample_birth(rng)};
		EXPECT_GT(b.e_impact, 0.0);
		EXPECT_GE(b.a_impact, 0.0);
		EXPECT_LE(b.a_impact, 90.0);
		const auto f {formula.sample_birth(rng)};
		EXPECT_EQ(f.a_impact, 60.0);
	}
}

// With a yield that only the fast ions get (a step at 1.2 times the
// formula's mean energy), the sheath model's yield is the flux-weighted
// share of its samples above the step, and the ions picked for births
// are all above it
TEST_F(SputteringTable, SheathImpactsStep)
{
	auto opts {slab_wall_options("ymin", 10.0)};
	opts.set_slab_flow(1e4);
	const auto bkg {Background::read_bkg(opts)};
	const Sputtering::Surface surface(bkg);
	const double e0 {1.2 * (Sputtering::sheath_energy(20.0, 20.0, 5.4858e-4,
		2.014, 1) + 40.0)};
	write(step_table(e0, 0.5));
	const Sputtering::YieldTable table(m_path);
	const Sputtering::Source sheath(surface, bkg, table, opts,
		Sheath::Model::brooks);
	ASSERT_NE(sheath.impacts(), nullptr);
	double above {};
	for (const auto& s : sheath.impacts()->samples(0))
		if (s.energy * 20.0 > e0) above += s.weight;
	ASSERT_GT(above, 0.0);
	EXPECT_NEAR(sheath.yield()[0], 0.5 * above, 1e-9);

	pcg32 rng(9, 1);
	for (int n {}; n < 200; ++n)
		EXPECT_GT(sheath.sample_birth(rng).e_impact, e0);
}

// An atom absorbed at the slab wall lands on the element under it
TEST(SlabWall, Deposition)
{
	const auto bkg {Background::read_bkg(slab_wall_options("ymin"))};
	const Sputtering::Surface surface(bkg);
	Sputtering::Deposition dep(surface, 74);
	const double side {0.05 / 8.0};
	const std::array<double, 3> pos {2.5 * side, -0.0251,
		-0.025 + 6.5 * side};
	dep.add_at_wall_face(bkg.get_zones(), 0, Zones::ymin, pos, 1e-6, 3,
		2.0);
	double total {};
	for (std::size_t e {}; e < surface.size(); ++e)
	{
		const auto& el {surface.elements()[e]};
		const double expected {el.i == 2 && el.k == 6 ? 2.0 : 0.0};
		EXPECT_DOUBLE_EQ(dep.redeposited()[e], expected);
		total += dep.redeposited()[e];
	}
	EXPECT_DOUBLE_EQ(total, 2.0);
	EXPECT_DOUBLE_EQ(dep.wall_face(), 0.0);
}

// What the slab won't do with a wall or a turned field
TEST(SlabWall, Checks)
{
	auto periodic {slab_wall_options("ymin")};
	periodic.set_max_ybound_type("periodic");
	EXPECT_EXIT(Background::read_bkg(periodic),
		testing::ExitedWithCode(EXIT_FAILURE), "can't be periodic");

	auto bent {slab_wall_options("ymin")};
	bent.set_slab_major_radius(1.0);
	EXPECT_EXIT(Background::read_bkg(bent),
		testing::ExitedWithCode(EXIT_FAILURE), "straight slab");

	auto e_field {slab_wall_options("ymin", 10.0)};
	e_field.set_slab_efield_y(100.0);
	EXPECT_EXIT(Background::read_bkg(e_field),
		testing::ExitedWithCode(EXIT_FAILURE), "slab_b_angle");

	auto too_far {slab_wall_options("ymin", 91.0)};
	EXPECT_EXIT(Background::read_bkg(too_far),
		testing::ExitedWithCode(EXIT_FAILURE), "-90 to 90");

	// Without a wall the slab has no surface to sputter
	EXPECT_EXIT(Sputtering::Surface {Background::read_bkg(
		slab_wall_options("none"))}, testing::ExitedWithCode(EXIT_FAILURE),
		"slab_wall");
}

// A return is prompt if it comes by prompt_by, and is tallied by element
// and charge as well as with the rest. One that never ionized (prompt_by
// below 0) or came later isn't.
TEST(SlabWall, PromptTally)
{
	const auto bkg {Background::read_bkg(slab_wall_options("ymin"))};
	const Sputtering::Surface surface(bkg);
	Sputtering::Deposition dep(surface, 74);
	const double side {0.05 / 8.0};
	const std::array<double, 3> pos {0.5 * side, -0.0251, -0.025
		+ 0.5 * side};
	const auto& zones {bkg.get_zones()};
	dep.add_at_wall_face(zones, 0, Zones::ymin, pos, 1e-6, 1, 1.0, 2e-6);
	dep.add_at_wall_face(zones, 0, Zones::ymin, pos, 2e-6, 2, 2.0, 2e-6);
	dep.add_at_wall_face(zones, 0, Zones::ymin, pos, 3e-6, 1, 4.0, 2e-6);
	dep.add_at_wall_face(zones, 0, Zones::ymin, pos, 1e-6, 0, 8.0, -1.0);
	dep.add_at_wall_face(zones, 0, Zones::ymin, pos, 1e-6, 1, 16.0);

	EXPECT_DOUBLE_EQ(dep.total(), 31.0);
	EXPECT_DOUBLE_EQ(dep.prompt_total(), 3.0);
	EXPECT_DOUBLE_EQ(dep.prompt_charge()[1], 1.0);
	EXPECT_DOUBLE_EQ(dep.prompt_charge()[2], 2.0);
	double on_elements {};
	for (std::size_t e {}; e < surface.size(); ++e)
	{
		const auto& el {surface.elements()[e]};
		const bool here {el.i == 0 && el.k == 0};
		EXPECT_DOUBLE_EQ(dep.prompt()[e], here ? 3.0 : 0.0);
		on_elements += dep.prompt()[e];
	}
	EXPECT_DOUBLE_EQ(on_elements, 3.0);
}

// The sheath field on a toroidal background, in front of the 45 degree
// staircase wall: B along Z is 45 degrees from the wall's normal
// (-1, -1) / sqrt(2) in (R, Z). In plasma cell (6, 3) the particle is 0.8
// of the way to the cell's xmax face, so 0.2 cell widths (2 mm) in R from
// it, 2 mm / sqrt(2) along the normal, nearer than the ymax face. The
// field is along the normal turned to the particle's toroidal angle.
TEST(SheathField, DiagonalWall)
{
	const auto bkg {diagonal_with_plasma(Plasma {})};
	const Sputtering::Surface surface(bkg);
	const Sheath::Sheath sheath(surface, bkg, Sheath::Model::brooks);
	ASSERT_EQ(sheath.size(), surface.size());

	const double phi {0.13};
	const std::array<double, 3> pos {6.8, 3.5, phi};
	std::size_t e {};
	double d {};
	ASSERT_TRUE(sheath.distance(0, 6, 3, 1, pos, e, d));
	const auto& el {surface.elements()[e]};
	EXPECT_EQ(el.face, Zones::xmax);
	EXPECT_EQ(el.i, 6);
	EXPECT_EQ(el.j, 3);
	EXPECT_NEAR(d, 0.002 / std::sqrt(2.0), 1e-12);

	const auto& p {sheath.params(0, e)};
	EXPECT_NEAR(p.theta, S3XTest::pi / 4.0, 1e-9);
	const double en {Sheath::field(p, Sheath::Model::brooks, d)};
	ASSERT_LT(en, 0.0);
	const auto E {sheath.field_at(0, 0, 6, 3, 1, pos)};
	const auto n {el.normal_xyz(phi, bkg.get_phi_sign())};
	for (int c {}; c < 3; ++c) EXPECT_NEAR(E[c], en * n[c], 1e-9 * -en);

	// Pulling toward the wall, up and out in R
	EXPECT_GT(E[0] * std::cos(phi) + E[1] * std::sin(phi), 0.0);
	EXPECT_GT(E[2], 0.0);

	// A cell away from the wall has no field
	const auto none {sheath.field_at(0, 0, 6, 1, 1, {6.5, 1.5, phi})};
	for (const double c : none) EXPECT_DOUBLE_EQ(c, 0.0);
}

// Impact tallies on the slab wall at y min (normal +Y): a particle coming
// straight in is at 0 degrees, one coming in at 45 degrees is at 45, and
// one at 60 degrees from the normal at 60. The energies are 0.5 m v^2,
// and the mean time since birth is weighted by rate.
TEST(SlabWall, ImpactTally)
{
	const auto bkg {Background::read_bkg(slab_wall_options("ymin"))};
	const Sputtering::Surface surface(bkg);
	Sputtering::Deposition dep(surface, 74);
	const double side {0.05 / 8.0};
	const std::array<double, 3> pos {0.5 * side, -0.0251, -0.025
		+ 0.5 * side};
	const auto& zones {bkg.get_zones()};
	const double mass {183.84 * 1.661e-27};
	const double v {1e4};
	const double energy {0.5 * mass * v * v / 1.602e-19};
	const double s {std::sqrt(0.5)};
	const Sputtering::Impact straight {{0.0, -v, 0.0}, mass};
	const Sputtering::Impact oblique {{v * s, -v * s, 0.0}, mass};
	const Sputtering::Impact steep {{0.0, -0.5 * v, v * std::sqrt(0.75)},
		mass};
	dep.add_at_wall_face(zones, 0, Zones::ymin, pos, 1e-6, 1, 1.0, -1.0,
		&straight);
	dep.add_at_wall_face(zones, 0, Zones::ymin, pos, 2e-6, 1, 2.0, -1.0,
		&oblique);
	dep.add_at_wall_face(zones, 0, Zones::ymin, pos, 4e-6, 1, 1.0, -1.0,
		&steep);

	// Without an Impact, only the time since birth counts
	dep.add_at_wall_face(zones, 0, Zones::ymin, pos, 0.0, 1, 4.0);

	EXPECT_NEAR(dep.mean_energy(), energy, 1e-9 * energy);
	EXPECT_NEAR(dep.mean_angle(), (0.0 + 2.0 * 45.0 + 60.0) / 4.0, 1e-9);
	EXPECT_NEAR(dep.mean_lifetime(), (1e-6 + 4e-6 + 4e-6) / 8.0, 1e-18);
	EXPECT_DOUBLE_EQ(dep.impact_angle()[0], 1.0);
	EXPECT_DOUBLE_EQ(dep.impact_angle()[44] + dep.impact_angle()[45], 2.0);
	EXPECT_DOUBLE_EQ(dep.impact_angle()[59] + dep.impact_angle()[60], 1.0);
	double in_bins {};
	for (const double r : dep.impact_energy()) in_bins += r;
	EXPECT_DOUBLE_EQ(in_bins, 4.0);
	const auto& edges {dep.energy_edges()};
	for (std::size_t b {}; b < dep.impact_energy().size(); ++b)
	{
		const bool here {edges[b] <= energy && energy < edges[b + 1]};
		EXPECT_DOUBLE_EQ(dep.impact_energy()[b], here ? 4.0 : 0.0);
	}
	double on_elements {};
	for (const double r : dep.element_energy()) on_elements += r;
	EXPECT_NEAR(on_elements, 4.0 * energy, 1e-9 * energy);
}

// Roulette keeps weights at or above the minimum, and below it keeps the
// expected weight: the particle carries on with the minimum or ends
TEST(SelfSputtering, Roulette)
{
	EXPECT_DOUBLE_EQ(Sputtering::roulette(2.0, 1.0, 0.99), 2.0);
	EXPECT_DOUBLE_EQ(Sputtering::roulette(1.0, 1.0, 0.99), 1.0);
	EXPECT_DOUBLE_EQ(Sputtering::roulette(0.0, 1.0, 0.0), 0.0);
	EXPECT_DOUBLE_EQ(Sputtering::roulette(0.3, 0.0, 0.99), 0.3);
	EXPECT_DOUBLE_EQ(Sputtering::roulette(0.3, 1.0, 0.29), 1.0);
	EXPECT_DOUBLE_EQ(Sputtering::roulette(0.3, 1.0, 0.31), 0.0);

	pcg32 rng(11, 3);
	const int n {1000000};
	const double w {0.3};
	double sum {};
	for (int k {}; k < n; ++k)
	{
		const double out {Sputtering::roulette(w, 1.0, rng.next_double())};
		ASSERT_TRUE(out == 0.0 || out == 1.0);
		sum += out;
	}
	// Binomial: standard error sqrt(w (1 - w) / n)
	EXPECT_NEAR(sum / n, w, 4.0 * std::sqrt(w * (1.0 - w) / n));
}

// The impact energy and angle: from the velocity alone for a neutral or
// with the sheath field on (phi_wall = 0), and with q phi_wall added to
// the motion along the normal for an ion without it
TEST(SelfSputtering, Impact)
{
	const double mass {183.84 * Constants::amu_to_kg};
	const double ev {Constants::ev_to_j};
	const double v {1e4};
	const double ke {0.5 * mass * v * v / ev};
	const std::array<double, 3> n {0.0, 1.0, 0.0};
	const double s {std::sqrt(0.5)};
	const std::array<double, 3> oblique {v * s, -v * s, 0.0};

	const auto neutral {Sputtering::self_impact(oblique, mass, n, 0, 50.0)};
	EXPECT_NEAR(neutral.energy, ke, 1e-9 * ke);
	EXPECT_NEAR(neutral.angle, 45.0, 1e-9);

	const auto on {Sputtering::self_impact(oblique, mass, n, 2, 0.0)};
	EXPECT_NEAR(on.energy, ke, 1e-9 * ke);
	EXPECT_NEAR(on.angle, 45.0, 1e-9);

	const auto off {Sputtering::self_impact(oblique, mass, n, 2, 50.0)};
	EXPECT_NEAR(off.energy, ke + 100.0, 1e-9 * (ke + 100.0));
	const double v_n {std::sqrt(0.5 * v * v + 2.0 * 2.0 * 50.0 * ev
		/ mass)};
	EXPECT_NEAR(off.angle, std::atan2(v * s, v_n) * 180.0 / S3XTest::pi,
		1e-9);

	// Moving away from the wall: its speed along the normal still counts
	const auto away {Sputtering::self_impact({v * s, v * s, 0.0}, mass, n,
		2, 50.0)};
	EXPECT_NEAR(away.energy, off.energy, 1e-9 * off.energy);
	EXPECT_NEAR(away.angle, off.angle, 1e-9);
}

namespace
{
	// A slab with a wall at y min, B at 10 degrees to it and a flow onto
	// it, its D+ source from a constant yield, and the self-sputtering
	// options
	Options::Options self_options()
	{
		auto opts {slab_wall_options("ymin", 10.0)};
		opts.set_slab_flow(1e4);
		opts.set_imp_num(100);
		opts.set_sput_self("on");
		opts.set_sput_self_min_weight(0.1);
		opts.set_sput_self_max_gen(3);
		return opts;
	}

	// A W on W table with a constant yield from 1 eV
	Table self_table(const double y)
	{
		Table t {constant_table(y)};
		t.projectile_Z = 74.0;
		t.projectile_mass_amu = 183.84;
		return t;
	}
}

// A W ion landing on the slab wall without the sheath field: its impact
// energy gets the element's floating drop, its weight becomes w Y, and it
// carries on as a neutral atom just inside the plasma cell where it hit,
// moving into the plasma no faster than Lambda E (Lambda = 1 for W on W)
// allows. The tallies get w Y.
TEST_F(SputteringTable, SelfRebirth)
{
	write(constant_table(0.02));
	write(self_table(0.5), m_self_path);
	const auto opts {self_options()};
	const auto bkg {Background::read_bkg(opts)};
	const Sputtering::Surface surface(bkg);
	const Sputtering::Source source(surface, bkg,
		Sputtering::YieldTable {m_path}, opts);
	Sputtering::SelfSource self(source, Sputtering::YieldTable {m_self_path},
		opts, false);
	const double w0 {source.source_rate() / 100.0};
	EXPECT_NEAR(self.min_weight(), 0.1 * w0, 1e-12 * w0);

	const double side {0.05 / 8.0};
	const std::array<double, 3> pos {2.3 * side, -0.0251, -0.025
		+ 5.6 * side};
	std::size_t e {surface.size()};
	for (std::size_t k {}; k < surface.size(); ++k)
	{
		const auto& el {surface.elements()[k]};
		if (el.i == 2 && el.k == 5) e = k;
	}
	ASSERT_LT(e, surface.size());

	const double mass {183.84 * Constants::amu_to_kg};
	const double ev {Constants::ev_to_j};
	const double v {1e4};
	const Sputtering::Impact impact {{0.0, -v, 0.0}, mass};
	const double drop {source.phi_wall()[e]};
	// The floating drop for D at Ti = Te: 2.84 Te
	EXPECT_NEAR(drop / 20.0, -0.5 * std::log(2.0 * S3XTest::pi
		* bkg.get_elec_mass_amu() / bkg.get_ion_mass_amu() * 2.0), 1e-9);
	const double energy {0.5 * mass * v * v / ev + 3.0 * drop};

	pcg32 rng(5, 2);
	const double y_cell {-0.025 + 1e-3 * side};
	for (int n {}; n < 200; ++n)
	{
		const auto r {self.rebirth(e, 0, pos, impact, 3, w0, 0, rng)};
		ASSERT_TRUE(r.alive);
		EXPECT_FALSE(r.avalanche);
		EXPECT_DOUBLE_EQ(r.weight, 0.5 * w0);
		EXPECT_EQ(r.zone, 0);
		EXPECT_NEAR(r.x, pos[0], 1e-15);
		EXPECT_NEAR(r.y, y_cell, 1e-12);
		EXPECT_NEAR(r.z, pos[2], 1e-15);
		EXPECT_EQ(r.i, 2);
		EXPECT_EQ(r.k, 5);
		EXPECT_GT(r.vY, 0.0);
		const double e_atom {0.5 * mass * (r.vX * r.vX + r.vY * r.vY
			+ r.vZ * r.vZ) / ev};
		EXPECT_LE(e_atom, energy * (1.0 + 1e-9));
	}
	EXPECT_NEAR(self.mean_energy(), energy, 1e-9 * energy);
	EXPECT_NEAR(self.mean_angle(), 0.0, 1e-9);
	EXPECT_NEAR(self.self_sputtered()[e], 200 * 0.5 * w0, 1e-9 * w0);
	EXPECT_NEAR(self.self_total(), 200 * 0.5 * w0, 1e-9 * w0);
	EXPECT_NEAR(self.generation_weight()[1], 200 * 0.5 * w0, 1e-9 * w0);
	EXPECT_NEAR(self.generation_landed()[0], 200 * w0, 1e-9 * w0);
	EXPECT_NEAR(self.mean_yield(), 0.5, 1e-12);

	// With the sheath field on, the velocity alone
	Sputtering::SelfSource with(source, Sputtering::YieldTable {m_self_path},
		opts, true);
	with.rebirth(e, 0, pos, impact, 3, w0, 0, rng);
	EXPECT_NEAR(with.mean_energy(), 0.5 * mass * v * v / ev, 1e-6);

	// A neutral below the table's lowest energy (1 eV) sputters nothing
	const Sputtering::Impact slow {{0.0, -100.0, 0.0}, mass};
	const auto none {with.rebirth(e, 0, pos, slow, 0, w0, 0, rng)};
	EXPECT_FALSE(none.alive);
	EXPECT_EQ(with.num_below_threshold(), 1);

	// Past sput_self_max_gen (3), the chain ends as an avalanche
	const auto last {self.rebirth(e, 0, pos, impact, 3, w0, 2, rng)};
	EXPECT_TRUE(last.alive);
	const auto over {self.rebirth(e, 0, pos, impact, 3, w0, 3, rng)};
	EXPECT_FALSE(over.alive);
	EXPECT_TRUE(over.avalanche);
}

// Light particles play roulette: they carry on with the minimum weight
// or end, keeping the expected weight, and what's lost is tallied
TEST_F(SputteringTable, SelfRoulette)
{
	write(constant_table(0.02));
	write(self_table(0.5), m_self_path);
	const auto opts {self_options()};
	const auto bkg {Background::read_bkg(opts)};
	const Sputtering::Surface surface(bkg);
	const Sputtering::Source source(surface, bkg,
		Sputtering::YieldTable {m_path}, opts);
	Sputtering::SelfSource self(source, Sputtering::YieldTable {m_self_path},
		opts, true);
	const double w_min {self.min_weight()};
	const double w {0.1 * w_min};	// w Y = 0.05 w_min
	const Sputtering::Impact impact {{0.0, -1e4, 0.0}, 183.84
		* Constants::amu_to_kg};
	const double side {0.05 / 8.0};
	pcg32 rng(9, 4);
	const int n {20000};
	int alive {};
	for (int k {}; k < n; ++k)
	{
		const auto r {self.rebirth(0, 0, {0.5 * side, -0.0251, -0.025
			+ 0.5 * side}, impact, 1, w, 0, rng)};
		if (r.alive)
		{
			EXPECT_DOUBLE_EQ(r.weight, w_min);
			++alive;
		}
	}
	const double p {0.05};
	EXPECT_NEAR(static_cast<double>(alive) / n, p, 4.0 * std::sqrt(p
		* (1.0 - p) / n));
	EXPECT_NEAR(self.generation_weight()[1], alive * w_min, 1e-9 * n * w);
	EXPECT_NEAR(self.roulette_lost(), (n - alive) * 0.5 * w, 1e-9 * n * w);
	EXPECT_NEAR(self.self_total(), n * 0.5 * w, 1e-9 * n * w);
}

// The self-sputtering table has to be the impurity on the D+ table's wall
TEST_F(SputteringTable, SelfChecks)
{
	write(constant_table(0.02));
	const auto opts {self_options()};
	const auto bkg {Background::read_bkg(opts)};
	const Sputtering::Surface surface(bkg);
	const Sputtering::Source source(surface, bkg,
		Sputtering::YieldTable {m_path}, opts);

	write(constant_table(0.5), m_self_path);
	EXPECT_EXIT(Sputtering::SelfSource(source,
		Sputtering::YieldTable {m_self_path}, opts, false),
		testing::ExitedWithCode(EXIT_FAILURE), "atomic number");

	Table light {self_table(0.5)};
	light.projectile_mass_amu = 150.0;
	write(light, m_self_path);
	EXPECT_EXIT(Sputtering::SelfSource(source,
		Sputtering::YieldTable {m_self_path}, opts, false),
		testing::ExitedWithCode(EXIT_FAILURE), "mass");
}

// sput_impact_opt: "formula" (the default) or "sheath", with its control
// integer; anything else stops the run
TEST(SputteringOptions, ImpactOpt)
{
	Options::Options opts {};
	EXPECT_EQ(opts.sput_impact_opt(), "formula");
	EXPECT_EQ(opts.sput_impact_opt_int(), 0);
	opts.set_sput_impact_opt("sheath");
	EXPECT_EQ(opts.sput_impact_opt(), "sheath");
	EXPECT_EQ(opts.sput_impact_opt_int(), 1);
	opts.set_sput_impact_opt("formula");
	EXPECT_EQ(opts.sput_impact_opt_int(), 0);
	EXPECT_EXIT(opts.set_sput_impact_opt("bogus"),
		testing::ExitedWithCode(EXIT_FAILURE), "sput_impact_opt");
}

// Through the boundary check: a W ion past the slab wall is redeposited
// and carries on as a neutral W atom of generation 1 in the plasma, with
// its time since birth reset
TEST_F(SputteringTable, SelfBoundary)
{
	write(constant_table(0.02));
	write(self_table(2.0), m_self_path);
	const auto opts {self_options()};
	const auto bkg {Background::read_bkg(opts)};
	const Sputtering::Surface surface(bkg);
	const Sputtering::Source source(surface, bkg,
		Sputtering::YieldTable {m_path}, opts);
	Sputtering::Deposition dep(surface, 74);
	Sputtering::SelfSource self(source, Sputtering::YieldTable {m_self_path},
		opts, false);
	dep.set_self(&self);

	Slots::Slots slots(1);
	slots.set_mass(183.84 * Constants::amu_to_kg);
	slots.set_Z(74);
	const double side {0.05 / 8.0};
	slots.set_x(0, 0.5 * side);
	slots.set_y(0, -0.0251);
	slots.set_z(0, -0.025 + 0.5 * side);
	slots.set_vY(0, -1e4);
	slots.set_q(0, 2);
	const double w0 {source.source_rate() / 100.0};
	slots.set_weight(0, w0);
	slots.set_pt(0, 3e-6);
	slots.set_t(0, 1e-6);
	slots.set_state(0, 0);

	std::vector<pcg32> rngs {};
	for (int k {}; k < omp_get_max_threads(); ++k) rngs.emplace_back(k + 1, 1);
	Boundary::Fates fates {};
	Boundary::check_bounds_cpu(slots, bkg, opts, rngs, fates, &dep);

	EXPECT_EQ(slots.state()[0], 0);
	EXPECT_EQ(slots.q()[0], 0);
	EXPECT_EQ(slots.generation()[0], 1);
	EXPECT_DOUBLE_EQ(slots.weight()[0], 2.0 * w0);
	EXPECT_DOUBLE_EQ(slots.pt()[0], 0.0);
	EXPECT_DOUBLE_EQ(slots.t()[0], 1e-6);
	EXPECT_GT(slots.y()[0], -0.025);
	EXPECT_GT(slots.vY()[0], 0.0);
	EXPECT_EQ(fates.wall, 1);
	EXPECT_EQ(fates.self_sputtered, 1);
	EXPECT_DOUBLE_EQ(dep.total(), w0);
	EXPECT_DOUBLE_EQ(self.self_total(), 2.0 * w0);
}
