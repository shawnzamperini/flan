/*
* Tests for reading a HERMES (JAX-Hermes3) background (read_hermes.cpp) on
* a small single-null grid and snapshots written by the test (see
* hermes_test_files.h): the six zones and their faces, the geometry and
* field, the unnormalized plasma, picking snapshots, and the errors.
*/

#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "background.h"
#include "constants.h"
#include "hermes_test_files.h"
#include "options.h"
#include "read_bkg.h"
#include "read_hermes.h"
#include "zones.h"

namespace
{
	constexpr int nz {4};
	constexpr int zperiod {2};
	constexpr double nnorm {2e19};
	constexpr double tnorm {40.0};
	constexpr double bnorm {1.5};
	constexpr double aa {2.0};

	// Global grid index of each zone's first cell (x0, y0) for the test
	// grid: ixseps1 = 5, jyseps1_1 = 1, jyseps2_2 = 5
	constexpr std::array<std::array<int, 2>, 6> origin {{
		{2, 0}, {2, 2}, {2, 6}, {5, 0}, {5, 2}, {5, 6}}};

	class ReadHermesTest : public testing::Test
	{
	protected:
		std::filesystem::path m_dir {std::filesystem::temp_directory_path()
			/ ("flan_test_hermes_" + std::to_string(getpid()))};
		HermesTest::TestGrid m_grid {};

		void SetUp() override
		{
			// Backgrounds use OpenMP, so death tests run in a new process
			GTEST_FLAG_SET(death_test_style, "threadsafe");
			std::filesystem::create_directories(m_dir / "snaps");
			write();
		}

		void TearDown() override {std::filesystem::remove_all(m_dir);}

		// Writes the grid and three snapshots at t = 100, 200, 300 (and a
		// file that isn't a snapshot)
		void write()
		{
			HermesTest::write_grid(grid_path(), m_grid);
			for (int t {}; t < 3; ++t)
			{
				HermesTest::write_snapshot((m_dir / "snaps").string(), m_grid,
					nz, 100 * (t + 1), t, 100.0 * (t + 1));
			}
			std::ofstream((m_dir / "snaps" / "notes.txt").string()) << "x";
		}

		std::string grid_path() const {return (m_dir / "grid.nc").string();}

		Options::Options options() const
		{
			Options::Options opts {};
			opts.set_bkg_source("hermes");
			opts.set_hermes_grid_file(grid_path());
			opts.set_hermes_snap_dir((m_dir / "snaps").string());
			opts.set_hermes_snap_start(0);
			opts.set_hermes_snap_end(-1);
			opts.set_hermes_nnorm(nnorm);
			opts.set_hermes_tnorm(tnorm);
			opts.set_hermes_bnorm(bnorm);
			opts.set_hermes_ion_mass_amu(aa);
			opts.set_hermes_zperiod(zperiod);
			return opts;
		}
	};
}

// Six zones with global grid indices as x and y, the toroidal wedge in z,
// and the single-null connections
TEST_F(ReadHermesTest, Zones)
{
	const auto bkg {Hermes::read_hermes(options())};
	const Zones::Zones& z {bkg.get_zones()};
	ASSERT_EQ(z.get_nzones(), 6);
	EXPECT_EQ(z.get_ghost(), 2);
	const std::array<std::array<int, 2>, 6> size {{
		{3, 2}, {3, 4}, {3, 2}, {3, 2}, {3, 4}, {3, 2}}};
	const double dz {2.0 * Constants::pi / (zperiod * nz)};
	for (int zone {}; zone < 6; ++zone)
	{
		ASSERT_EQ(z.get_nx(zone), size[zone][0]);
		ASSERT_EQ(z.get_ny(zone), size[zone][1]);
		ASSERT_EQ(z.get_nz(zone), nz);
		EXPECT_EQ(z.get_grid_x(zone).front(), origin[zone][0]);
		EXPECT_EQ(z.get_grid_y(zone).front(), origin[zone][1]);
		for (int k {}; k < nz; ++k)
			EXPECT_NEAR(z.get_z(zone)[k], k * dz, 1e-14);
		EXPECT_EQ(z.get_face(zone, Zones::zmin), zone);
		EXPECT_EQ(z.get_face(zone, Zones::zmax), zone);
	}

	using namespace Zones;
	auto faces = [&](int zone) {return std::array<int, 4> {
		z.get_face(zone, xmin), z.get_face(zone, xmax),
		z.get_face(zone, ymin), z.get_face(zone, ymax)};};
	EXPECT_EQ(faces(0), (std::array<int, 4> {face_absorbing, 3, face_wall,
		2}));
	EXPECT_EQ(faces(1), (std::array<int, 4> {face_core, 4, 1, 1}));
	EXPECT_EQ(faces(2), (std::array<int, 4> {face_absorbing, 5, 0,
		face_wall}));
	EXPECT_EQ(faces(3), (std::array<int, 4> {0, face_absorbing, face_wall,
		4}));
	EXPECT_EQ(faces(4), (std::array<int, 4> {1, face_absorbing, 3, 5}));
	EXPECT_EQ(faces(5), (std::array<int, 4> {2, face_absorbing, 4,
		face_wall}));
}

// Cell centers from the grid's corners in a right-handed (R, z, Z), and
// the field from Brxy, Bzxy and Btxy
TEST_F(ReadHermesTest, GeometryAndField)
{
	const auto bkg {Hermes::read_hermes(options())};
	EXPECT_EQ(bkg.get_phi_sign(), 1.0);
	const Zones::Zones& z {bkg.get_zones()};
	const auto& g {m_grid};
	for (int zone {}; zone < 6; ++zone)
	for (int i {}; i < z.get_nx(zone); ++i)
	for (int j {}; j < z.get_ny(zone); ++j)
	for (int k {}; k < nz; ++k)
	{
		const double phi {z.get_z(zone)[k]};
		const double R {g.R0 + (origin[zone][0] + i + 0.5) * g.dR};
		const double Z {g.Z0 + (origin[zone][1] + j + 0.5) * g.dZ};
		EXPECT_NEAR(bkg.get_X()(zone, i, j, k), R * std::cos(phi), 1e-12);
		EXPECT_NEAR(bkg.get_Y()(zone, i, j, k), R * std::sin(phi), 1e-12);
		EXPECT_NEAR(bkg.get_Z()(zone, i, j, k), Z, 1e-12);

		// B = Br R_hat + Bt phi_hat + Bz Z_hat
		const double bX {bkg.get_bX()(zone, 0, i, j, k)};
		const double bY {bkg.get_bY()(zone, 0, i, j, k)};
		EXPECT_NEAR(bX * std::cos(phi) + bY * std::sin(phi), g.Br, 1e-12);
		EXPECT_NEAR(-bX * std::sin(phi) + bY * std::cos(phi), g.Bt, 1e-12);
		EXPECT_NEAR(bkg.get_bZ()(zone, 0, i, j, k), g.Bz, 1e-12);
	}
}

// Every cell gets the unnormalized plasma of its global grid cell, and the
// times are t m_p / (e Bnorm)
TEST_F(ReadHermesTest, Plasma)
{
	const auto bkg {Hermes::read_hermes(options())};
	const Zones::Zones& z {bkg.get_zones()};
	const double cs0 {std::sqrt(Hermes::qe * tnorm / Hermes::mp)};
	const double t_unit {Hermes::mp / (Hermes::qe * bnorm)};
	const auto& times {bkg.get_times()};
	ASSERT_EQ(times.size(), 3u);
	for (int t {}; t < 3; ++t)
		EXPECT_NEAR(times[t], 100.0 * (t + 1) * t_unit, 1e-20);
	EXPECT_EQ(bkg.get_ion_mass_amu(), aa);

	using HermesTest::field_value;
	for (int zone {}; zone < 6; ++zone)
	for (int t {}; t < 3; ++t)
	for (int i {}; i < z.get_nx(zone); ++i)
	for (int j {}; j < z.get_ny(zone); ++j)
	for (int k {}; k < nz; ++k)
	{
		const int gi {origin[zone][0] + i};
		const int gj {origin[zone][1] + j};
		const double Ne {field_value(0, t, gi, gj, k)};
		auto rel = [](double a, double b) {return std::abs(a / b - 1.0);};
		EXPECT_LT(rel(bkg.get_ne()(zone, t, i, j, k), Ne * nnorm), 1e-14);
		EXPECT_LT(rel(bkg.get_te()(zone, t, i, j, k),
			field_value(1, t, gi, gj, k) / Ne * tnorm), 1e-14);
		EXPECT_LT(rel(bkg.get_ti()(zone, t, i, j, k),
			field_value(2, t, gi, gj, k) / Ne * tnorm), 1e-14);
		EXPECT_LT(rel(bkg.get_vp()(zone, t, i, j, k),
			field_value(4, t, gi, gj, k) * tnorm), 1e-14);

		// The flow along b is the parallel flow, since the drifts are
		// perpendicular to b
		const std::array<double, 3> b {bkg.get_bX()(zone, t, i, j, k),
			bkg.get_bY()(zone, t, i, j, k), bkg.get_bZ()(zone, t, i, j, k)};
		const double bmag {std::sqrt(b[0] * b[0] + b[1] * b[1]
			+ b[2] * b[2])};
		const double u_par {(bkg.get_uX()(zone, t, i, j, k) * b[0]
			+ bkg.get_uY()(zone, t, i, j, k) * b[1]
			+ bkg.get_uZ()(zone, t, i, j, k) * b[2]) / bmag};
		const double want {field_value(3, t, gi, gj, k) / (aa * Ne) * cs0};
		EXPECT_NEAR(u_par, want, 1e-9 * std::abs(want));
	}
}

// bkg_source = "hermes" goes through read_bkg
TEST_F(ReadHermesTest, ReadBkg)
{
	const auto bkg {Background::read_bkg(options())};
	EXPECT_EQ(bkg.get_nzones(), 6);
	EXPECT_EQ(bkg.get_times().size(), 3u);
}

// With the poloidal field against +y, the parallel flow along B flips
// sign; with a mix the direction is ambiguous and the run stops
TEST_F(ReadHermesTest, ParallelSign)
{
	EXPECT_EQ(Hermes::parallel_sign(Hermes::read_grid(grid_path())), 1.0);
	m_grid.Bz = -m_grid.Bz;
	HermesTest::write_grid(grid_path(), m_grid);
	EXPECT_EQ(Hermes::parallel_sign(Hermes::read_grid(grid_path())), -1.0);
	m_grid.bz_flip_row = 3;
	HermesTest::write_grid(grid_path(), m_grid);
	EXPECT_EXIT(Hermes::parallel_sign(Hermes::read_grid(grid_path())),
		testing::ExitedWithCode(EXIT_FAILURE), "ambiguous");
}

// Snapshots are picked by index in name order, with a stride and an end
// counted from the back (and a start, too)
TEST_F(ReadHermesTest, SelectSnapshots)
{
	const std::string dir {(m_dir / "snaps").string()};
	auto names = [&](int start, int end, int stride)
	{
		std::vector<std::string> out {};
		for (const auto& p : Hermes::select_snapshots(dir, start, end,
			stride)) out.push_back(std::filesystem::path(p).filename());
		return out;
	};
	EXPECT_EQ(names(0, -1, 1), (std::vector<std::string> {
		"snap0000100.npz", "snap0000200.npz", "snap0000300.npz"}));
	EXPECT_EQ(names(0, 2, 2), (std::vector<std::string> {
		"snap0000100.npz", "snap0000300.npz"}));
	EXPECT_EQ(names(1, -1, 1), (std::vector<std::string> {
		"snap0000200.npz", "snap0000300.npz"}));
	EXPECT_EQ(names(-2, -1, 1), (std::vector<std::string> {
		"snap0000200.npz", "snap0000300.npz"}));
	EXPECT_EXIT(Hermes::select_snapshots(dir, -4, -1, 1),
		testing::ExitedWithCode(EXIT_FAILURE), "has 3 snapshots");

	EXPECT_EXIT(Hermes::select_snapshots(dir, 0, 1, 2),
		testing::ExitedWithCode(EXIT_FAILURE), "at least two");
	EXPECT_EXIT(Hermes::select_snapshots(dir, 0, 3, 1),
		testing::ExitedWithCode(EXIT_FAILURE), "has 3 snapshots");
	EXPECT_EXIT(Hermes::select_snapshots(dir, 0, 2, 0),
		testing::ExitedWithCode(EXIT_FAILURE), "at least 1");
}

// Things the reader can't use stop the run with a message saying why
TEST_F(ReadHermesTest, Errors)
{
	auto opts {options()};
	opts.set_hermes_format("bout");
	EXPECT_EXIT(Hermes::read_hermes(opts),
		testing::ExitedWithCode(EXIT_FAILURE), "isn't implemented yet");

	opts = options();
	opts.set_hermes_grid_file("undefined");
	EXPECT_EXIT(Hermes::read_hermes(opts),
		testing::ExitedWithCode(EXIT_FAILURE), "must be set");

	opts = options();
	opts.set_hermes_zperiod(0);
	EXPECT_EXIT(Hermes::read_hermes(opts),
		testing::ExitedWithCode(EXIT_FAILURE), "at least 1");

	// A double null grid
	auto dn {m_grid};
	dn.jyseps1_2 = 5;
	HermesTest::write_grid(grid_path(), dn);
	EXPECT_EXIT(Hermes::read_hermes(options()),
		testing::ExitedWithCode(EXIT_FAILURE), "isn't a single null");
	HermesTest::write_grid(grid_path(), m_grid);

	// A snapshot with a different nz from the first
	HermesTest::write_snapshot((m_dir / "snaps").string(), m_grid, nz + 1,
		400, 3, 400.0);
	EXPECT_EXIT(Hermes::read_hermes(options()),
		testing::ExitedWithCode(EXIT_FAILURE), "first snapshot's nz");
	std::filesystem::remove(m_dir / "snaps" / "snap0000400.npz");

	// Times that don't increase
	HermesTest::write_snapshot((m_dir / "snaps").string(), m_grid, nz,
		400, 3, 50.0);
	EXPECT_EXIT(Hermes::read_hermes(options()),
		testing::ExitedWithCode(EXIT_FAILURE), "times don't increase");
}
