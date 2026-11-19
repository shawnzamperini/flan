/*
* Tests for reading a HERMES background (read_hermes.cpp,
* read_hermes_bout.cpp) on a small single-null grid with JAX-Hermes3
* snapshots and BOUT++ Hermes-3 dump files written by the test (see
* hermes_test_files.h): the six zones and their faces, the geometry and
* field, the unnormalized plasma, picking snapshots, the BOUT++ reader
* against the JAX one on the same plasma, the sputtering surface at the
* targets, and the errors.
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <tuple>
#include <unistd.h>
#include <vector>

#include "background.h"
#include "constants.h"
#include "hermes_test_files.h"
#include "options.h"
#include "read_bkg.h"
#include "read_hermes.h"
#include "sputtering.h"
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

		std::string run_dir() const {return (m_dir / "run").string();}

		// Writes BOUT++ dump files of the snapshots' plasma in run_dir,
		// replacing any there
		void write_dumps(const HermesTest::DumpSpec& spec = {}) const
		{
			std::filesystem::remove_all(run_dir());
			std::filesystem::create_directories(run_dir());
			HermesTest::write_bout_dumps(run_dir(), m_grid, spec);
		}

		// Options for reading the dumps. The normalization and zperiod
		// options are left at their defaults, which differ from the
		// dumps', since they must come from the dumps. The inductive
		// E_par is off, so the background matches the JAX-Hermes3 one
		// (which has no Apar); BoutInductiveEpar tests it.
		Options::Options bout_options() const
		{
			Options::Options opts {};
			opts.set_bkg_source("hermes");
			opts.set_hermes_format("bout");
			opts.set_hermes_grid_file(grid_path());
			opts.set_hermes_snap_dir(run_dir());
			opts.set_hermes_snap_start(0);
			opts.set_hermes_snap_end(-1);
			opts.set_hermes_ion_mass_amu(aa);
			opts.set_hermes_inductive_epar("off");
			return opts;
		}
	};

	// Largest |a / b - 1| over two fields (0 where both are 0), and
	// whether they have the same size
	double max_rel_diff(const Vectors::ZonedVector<BkgFPType>& a,
		const Vectors::ZonedVector<BkgFPType>& b)
	{
		const auto& x {a.get_data()};
		const auto& y {b.get_data()};
		if (x.size() != y.size()) return INFINITY;
		double m {};
		for (std::size_t i {}; i < x.size(); ++i)
		{
			if (x[i] == y[i]) continue;
			const double d {std::abs(x[i] - y[i])
				/ std::max(std::abs(x[i]), std::abs(y[i]))};
			m = std::isnan(d) ? INFINITY : std::max(m, d);
		}
		return m;
	}

	// Checks that two backgrounds have the same plasma, fields and flows:
	// exactly, except the flows (and what depends on them) to within tol
	void expect_same(const Background::Background& a,
		const Background::Background& b, const double tol)
	{
		EXPECT_EQ(a.get_times(), b.get_times());
		EXPECT_EQ(a.get_ion_mass_amu(), b.get_ion_mass_amu());
		EXPECT_EQ(a.get_zones().get_nz(0), b.get_zones().get_nz(0));
		EXPECT_EQ(max_rel_diff(a.get_ne(), b.get_ne()), 0.0);
		EXPECT_EQ(max_rel_diff(a.get_te(), b.get_te()), 0.0);
		EXPECT_EQ(max_rel_diff(a.get_ti(), b.get_ti()), 0.0);
		EXPECT_EQ(max_rel_diff(a.get_vp(), b.get_vp()), 0.0);
		EXPECT_EQ(max_rel_diff(a.get_eX(), b.get_eX()), 0.0);
		EXPECT_EQ(max_rel_diff(a.get_eY(), b.get_eY()), 0.0);
		EXPECT_EQ(max_rel_diff(a.get_eZ(), b.get_eZ()), 0.0);
		EXPECT_EQ(max_rel_diff(a.get_X(), b.get_X()), 0.0);
		EXPECT_EQ(max_rel_diff(a.get_bmag(), b.get_bmag()), 0.0);
		EXPECT_EQ(max_rel_diff(a.get_q_par_e(), b.get_q_par_e()), 0.0);
		EXPECT_EQ(max_rel_diff(a.get_q_par_i(), b.get_q_par_i()), 0.0);
		EXPECT_LE(max_rel_diff(a.get_uX(), b.get_uX()), tol);
		EXPECT_LE(max_rel_diff(a.get_uY(), b.get_uY()), tol);
		EXPECT_LE(max_rel_diff(a.get_uZ(), b.get_uZ()), tol);
	}
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

// The sputtering surface of a HERMES background is its four target faces
// (the ymin faces of the inner PFR and leg, the ymax faces of the outer
// ones), one element per cell along them and toroidal cell. The test
// grid's cells are rectangles in (R, Z), so the normals point along Z into
// the plasma and each element's area is dR x R x dphi.
TEST_F(ReadHermesTest, SputteringTargets)
{
	const auto bkg {Hermes::read_hermes(options())};
	const Sputtering::Surface surface(bkg);
	ASSERT_TRUE(surface.from_wall_faces());

	const std::array<std::array<int, 2>, 4> targets {{
		{0, Zones::ymin}, {2, Zones::ymax}, {3, Zones::ymin},
		{5, Zones::ymax}}};
	ASSERT_EQ(surface.size(), static_cast<std::size_t>(4 * 3 * nz));

	const double dphi {2.0 * Constants::pi / (zperiod * nz)};
	const HermesTest::TestGrid& g {m_grid};
	double total {};
	for (const auto& e : surface.elements())
	{
		const bool listed {std::find(targets.begin(), targets.end(),
			std::array<int, 2> {e.zone, static_cast<int>(e.face)})
			!= targets.end()};
		EXPECT_TRUE(listed) << "zone " << e.zone << " face " << e.face;
		const bool lower {e.face == Zones::ymin};
		EXPECT_EQ(e.j, lower ? 0 : bkg.get_zones().get_ny(e.zone) - 1);

		// Global x index of the cell, and the target's Z: the bottom of
		// global row 0 or the top of the last row
		const int gi {origin[e.zone][0] + e.i};
		const double R_mid {g.R0 + (gi + 0.5) * g.dR};
		EXPECT_NEAR(e.R_mid, R_mid, 1e-12);
		EXPECT_NEAR(e.Z_mid, lower ? g.Z0 : g.Z0 + g.ny * g.dZ, 1e-12);
		EXPECT_NEAR(e.n_R, 0.0, 1e-12);
		EXPECT_NEAR(e.n_Z, lower ? 1.0 : -1.0, 1e-12);
		EXPECT_NEAR(e.area, g.dR * R_mid * dphi, 1e-12);
		total += e.area;
	}
	EXPECT_NEAR(surface.total_area(), total, 1e-12);
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
		testing::ExitedWithCode(EXIT_FAILURE), "has no BOUT.dmp.0.nc");

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

// ---------------------------------------------------------------------------
// BOUT++ Hermes-3 dump files
// ---------------------------------------------------------------------------

// The same plasma written as BOUT++ dumps (split over 2 x 4 processes,
// NaN in the guard cells the reader must skip, an electromagnetic run with
// V written) and as JAX snapshots gives the same background, with the
// normalization and zperiod taken from the dumps
TEST_F(ReadHermesTest, BoutMatchesJax)
{
	const Options::Options defaults {};
	ASSERT_NE(defaults.hermes_nnorm(), nnorm);
	ASSERT_NE(defaults.hermes_tnorm(), tnorm);
	ASSERT_NE(defaults.hermes_bnorm(), bnorm);
	ASSERT_NE(defaults.hermes_zperiod(), zperiod);

	write_dumps();
	const auto jax {Hermes::read_hermes(options())};
	const auto bout {Hermes::read_hermes(bout_options())};
	expect_same(bout, jax, 0.0);
}

// Without V in the dumps the flow comes from NV: (NV - N Apar) / (AA N) in
// an electromagnetic run (NV is the canonical momentum there), NV / (AA N)
// without Apar. Rounding in the subtraction allows a tiny difference.
TEST_F(ReadHermesTest, BoutFlowFromNV)
{
	const auto jax {Hermes::read_hermes(options())};
	for (const bool apar : {true, false})
	{
		SCOPED_TRACE(apar ? "with Apar" : "without Apar");
		HermesTest::DumpSpec spec {};
		spec.write_v = false;
		spec.write_apar = apar;
		write_dumps(spec);
		expect_same(Hermes::read_hermes(bout_options()), jax, 1e-13);
	}
}

// Output times are picked like snapshots: by index, with a stride and
// negative indices counting from the end
TEST_F(ReadHermesTest, BoutTimeSelection)
{
	HermesTest::DumpSpec spec {};
	spec.ntimes = 5;
	write_dumps(spec);
	const double t_unit {Hermes::mp / (Hermes::qe * bnorm)};
	auto check = [&](int start, int end, int stride,
		const std::vector<int>& want)
	{
		SCOPED_TRACE(std::to_string(start) + ", " + std::to_string(end)
			+ ", " + std::to_string(stride));
		auto opts {bout_options()};
		opts.set_hermes_snap_start(start);
		opts.set_hermes_snap_end(end);
		opts.set_hermes_snap_stride(stride);
		const auto bkg {Hermes::read_hermes(opts)};
		const auto& times {bkg.get_times()};
		ASSERT_EQ(times.size(), want.size());
		for (std::size_t t {}; t < want.size(); ++t)
		{
			EXPECT_NEAR(times[t], 100.0 * (want[t] + 1) * t_unit, 1e-20);
			// Main SOL (zone 4), its first cell
			const double Ne {HermesTest::field_value(0, want[t],
				origin[4][0], origin[4][1], 0)};
			EXPECT_DOUBLE_EQ(bkg.get_ne()(4, t, 0, 0, 0), Ne * nnorm);
		}
	};
	check(0, -1, 2, {0, 2, 4});
	check(-2, -1, 1, {3, 4});
	check(1, 3, 1, {1, 2, 3});
}

// The normalization options aren't used for BOUT++ dumps, even if set
TEST_F(ReadHermesTest, BoutIgnoresNormOptions)
{
	write_dumps();
	auto opts {bout_options()};
	opts.set_hermes_nnorm(5e19);
	opts.set_hermes_tnorm(10.0);
	opts.set_hermes_bnorm(2.5);
	opts.set_hermes_zperiod(7);
	expect_same(Hermes::read_hermes(opts),
		Hermes::read_hermes(bout_options()), 0.0);
}

// Dumps the reader can't use stop the run with a message saying why
// hermes_inductive_epar = "on" adds -dApar/dt along B to E, with
// dApar/dt the centered difference over the output times on either side
// of each loaded one (one sided at the first and last), whatever the
// stride. E x B, and so the flow, doesn't change.
TEST_F(ReadHermesTest, BoutInductiveEpar)
{
	HermesTest::DumpSpec spec {};
	spec.ntimes = 6;
	write_dumps(spec);
	const double cs0 {std::sqrt(Hermes::qe * tnorm / Hermes::mp)};
	const double sign {Hermes::parallel_sign(Hermes::read_grid(
		grid_path()))};
	for (const auto& [start, stride, want] : {
		std::tuple {0, 2, std::vector<int> {0, 2, 4}},
		std::tuple {1, 2, std::vector<int> {1, 3, 5}},
		std::tuple {0, 1, std::vector<int> {0, 1, 2, 3, 4, 5}}})
	{
		SCOPED_TRACE("start " + std::to_string(start) + ", stride "
			+ std::to_string(stride));
		auto opts {bout_options()};
		opts.set_hermes_snap_start(start);
		opts.set_hermes_snap_stride(stride);
		const auto off {Hermes::read_hermes(opts)};
		opts.set_hermes_inductive_epar("on");
		const auto on {Hermes::read_hermes(opts)};
		EXPECT_EQ(max_rel_diff(on.get_uX(), off.get_uX()), 0.0);
		EXPECT_EQ(max_rel_diff(on.get_uZ(), off.get_uZ()), 0.0);

		const Zones::Zones& z {on.get_zones()};
		using Vec3 = std::array<double, 3>;
		auto dot = [](const Vec3& a, const Vec3& b)
			{return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];};
		for (int zone {}; zone < 6; ++zone)
		for (std::size_t t {}; t < want.size(); ++t)
		for (int i {}; i < z.get_nx(zone); ++i)
		for (int j {}; j < z.get_ny(zone); ++j)
		for (int k {}; k < z.get_nz(zone); ++k)
		{
			const int gi {origin[zone][0] + i};
			const int gj {origin[zone][1] + j};
			const int lo {std::max(want[t] - 1, 0)};
			const int hi {std::min(want[t] + 1, spec.ntimes - 1)};
			const double dadt {(HermesTest::apar_value(hi, gi, gj, k)
				- HermesTest::apar_value(lo, gi, gj, k)) / (100.0 * (hi
				- lo))};
			const double expected {-sign * dadt * bnorm * cs0};

			const Vec3 B {on.get_bX()(zone, t, i, j, k),
				on.get_bY()(zone, t, i, j, k), on.get_bZ()(zone, t, i, j, k)};
			const double Bm {std::sqrt(dot(B, B))};
			const Vec3 dE {on.get_eX()(zone, t, i, j, k)
				- off.get_eX()(zone, t, i, j, k),
				on.get_eY()(zone, t, i, j, k) - off.get_eY()(zone, t, i, j,
				k), on.get_eZ()(zone, t, i, j, k) - off.get_eZ()(zone, t, i,
				j, k)};
			const double par {dot(dE, B) / Bm};
			ASSERT_NEAR(par, expected, 1e-9 * std::abs(expected));
			for (int c {}; c < 3; ++c)
				ASSERT_NEAR(dE[c], par * B[c] / Bm, 1e-9 * std::abs(par));
		}
	}
}

// The inductive E_par is on by default, and a run without Apar (an
// electrostatic BOUT++ run, or JAX-Hermes3) goes on without it
TEST_F(ReadHermesTest, InductiveEparDefault)
{
	EXPECT_EQ(Options::Options {}.hermes_inductive_epar(), "on");

	auto jax {options()};
	jax.set_hermes_inductive_epar("on");
	const auto jax_on {Hermes::read_hermes(jax)};
	jax.set_hermes_inductive_epar("off");
	expect_same(jax_on, Hermes::read_hermes(jax), 0.0);

	HermesTest::DumpSpec spec {};
	spec.write_apar = false;
	write_dumps(spec);
	auto bout {bout_options()};
	const auto bout_off {Hermes::read_hermes(bout)};
	bout.set_hermes_inductive_epar("on");
	expect_same(Hermes::read_hermes(bout), bout_off, 0.0);
}

TEST_F(ReadHermesTest, BoutErrors)
{
	auto expect_error = [&](const HermesTest::DumpSpec& spec,
		const std::string& message, const Options::Options& opts)
	{
		SCOPED_TRACE(message);
		write_dumps(spec);
		EXPECT_EXIT(Hermes::read_hermes(opts),
			testing::ExitedWithCode(EXIT_FAILURE), message);
	};
	const auto opts {bout_options()};
	HermesTest::DumpSpec spec {};

	// A missing file
	write_dumps();
	std::filesystem::remove(std::filesystem::path(run_dir())
		/ "BOUT.dmp.3.nc");
	EXPECT_EXIT(Hermes::read_hermes(opts),
		testing::ExitedWithCode(EXIT_FAILURE), "are missing");

	spec = {};
	spec.nzpe = 2;
	expect_error(spec, "splits z", opts);

	spec = {};
	spec.bad_times_file = 2;
	expect_error(spec, "different output times", opts);

	spec = {};
	spec.ixseps1_shift = 1;
	expect_error(spec, "aren't on the grid", opts);

	spec = {};
	spec.duplicate_place = true;
	expect_error(spec, "unfilled", opts);

	spec = {};
	spec.ion = "d+";
	expect_error(spec, "hermes_ion_species", opts);

	spec = {};
	spec.write_v = false;
	spec.write_nv = false;
	expect_error(spec, "neither", opts);

	auto late {bout_options()};
	late.set_hermes_snap_end(3);
	expect_error({}, "has 3 output times", late);
}

// With hermes_ion_species, the ion variables can have another name
TEST_F(ReadHermesTest, BoutIonSpecies)
{
	HermesTest::DumpSpec spec {};
	spec.ion = "d+";
	write_dumps(spec);
	auto opts {bout_options()};
	opts.set_hermes_ion_species("d+");
	expect_same(Hermes::read_hermes(opts), Hermes::read_hermes(options()),
		0.0);
}
