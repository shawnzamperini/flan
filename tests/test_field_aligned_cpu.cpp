/*
* Tests for field-aligned interpolation and derivatives on a
* shifted-metric (HERMES) background: the zShift tables and their ghost
* cells, agreement of zShift with the field and geometry, the spectral
* shift, interpolation of scalars and vectors along the field line, and
* the field-aligned y derivative and E.
*
* The background is the small single-null test grid of
* hermes_test_files.h, with a strong toroidal field so the field line
* moves 3/4 of the wedge (a quarter of the torus) per y cell and three
* times around the wedge going around the core. The plasma is a function
* of zeta = z - zShift, i.e., constant along the field lines. Such a
* function only closes on itself around the core if the core's twist is
* a whole number of periods of the function, which real turbulence
* satisfies by itself; here the twist is made the same in every column.
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <gtest/gtest.h>
#include <numbers>
#include <string>
#include <unistd.h>
#include <vector>

#include "axisym_mesh.h"
#include "background.h"
#include "braginskii.h"
#include "hermes_test_files.h"
#include "options.h"
#include "read_hermes.h"
#include "utilities.h"
#include "zones.h"

namespace
{
	constexpr int nz {32};
	constexpr int zperiod {4};
	constexpr double period {2.0 * std::numbers::pi / zperiod};
	constexpr double tnorm {20.0};

	// Global grid index of each zone's first cell, as in read_hermes
	// (ixseps1 = 5, jyseps1_1 = 1, jyseps2_2 = 5)
	constexpr std::array<std::array<int, 2>, 6> origin {{
		{2, 0}, {2, 2}, {2, 6}, {5, 0}, {5, 2}, {5, 6}}};

	// The test grid with no radial field (so b . grad f has no x part)
	// and a twist of -3/4 of the wedge per y cell, so -3 wedges over the
	// 4 core rows
	HermesTest::TestGrid aligned_grid()
	{
		HermesTest::TestGrid g {};
		g.Br = 0.0;
		g.twist = -0.75 * period;
		return g;
	}

	// A field-aligned profile, periodic in the wedge
	double profile(double zeta)
	{
		return std::cos(2.0 * std::numbers::pi * zeta / period);
	}

	class FieldAlignedTest : public testing::Test
	{
	protected:
		std::filesystem::path m_dir {std::filesystem::temp_directory_path()
			/ ("flan_test_aligned_" + std::to_string(getpid()))};
		HermesTest::TestGrid m_grid {aligned_grid()};

		void SetUp() override
		{
			GTEST_FLAG_SET(death_test_style, "threadsafe");
			std::filesystem::create_directories(m_dir / "snaps");
			HermesTest::write_grid((m_dir / "grid.nc").string(), m_grid);

			// n = 1 + 0.5 f(zeta), T = Tnorm, no flow and
			// phi = 0.3 Tnorm f(zeta), at two times
			const double dz {period / nz};
			const auto& g {m_grid};
			auto fn = [&g, dz](int f, int i, int j, int k)
			{
				const double p {profile(k * dz - g.z_shift(i, j))};
				const double ne {1.0 + 0.5 * p};
				switch (f)
				{
					case 0: case 1: case 2: return ne;
					case 3: return 0.0;
					default: return 0.3 * p;
				}
			};
			for (int t {}; t < 2; ++t)
			{
				HermesTest::write_snapshot((m_dir / "snaps").string(),
					m_grid, nz, t, t, 10.0 * t, fn);
			}
		}

		void TearDown() override {std::filesystem::remove_all(m_dir);}

		Background::Background read() const
		{
			Options::Options opts {};
			opts.set_bkg_source("hermes");
			opts.set_hermes_grid_file((m_dir / "grid.nc").string());
			opts.set_hermes_snap_dir((m_dir / "snaps").string());
			opts.set_hermes_snap_start(0);
			opts.set_hermes_snap_end(-1);
			opts.set_hermes_tnorm(tnorm);
			opts.set_hermes_zperiod(zperiod);
			return Hermes::read_hermes(opts);
		}
	};
}

// Every cell has the grid's zShift. Ghost cells across a y face continue
// the field line from the zone's edge, including the core's twist across
// the branch cut and the jump from one PFR leg to the other. Ghost cells
// across an x face have the grid's zShift of their cell.
TEST_F(FieldAlignedTest, ZShiftTables)
{
	const auto bkg {read()};
	const Zones::Zones& z {bkg.get_zones()};
	const auto& g {m_grid};
	for (int zone {}; zone < 6; ++zone)
	{
		ASSERT_TRUE(z.has_z_shift(zone));
		const int nx {z.get_nx(zone)};
		const int ny {z.get_ny(zone)};
		for (int i {}; i < nx; ++i)
		{
			const int gi {origin[zone][0] + i};
			const double s {g.twist_per_cell(gi)};
			for (int j {}; j < ny; ++j)
			{
				EXPECT_DOUBLE_EQ(z.z_shift(zone, i, j),
					g.z_shift(gi, origin[zone][1] + j));
			}
			for (int d : {1, 2})
			{
				if (z.get_face(zone, Zones::ymin) >= 0)
				{
					EXPECT_NEAR(z.z_shift(zone, i, -d), z.z_shift(zone, i, 0)
						- d * s, 1e-12) << "zone " << zone << " i " << i;
				}
				if (z.get_face(zone, Zones::ymax) >= 0)
				{
					EXPECT_NEAR(z.z_shift(zone, i, ny - 1 + d),
						z.z_shift(zone, i, ny - 1) + d * s, 1e-12)
						<< "zone " << zone << " i " << i;
				}
			}
		}
		for (int i : {-2, -1, nx, nx + 1})
		{
			const auto face {i < 0 ? Zones::xmin : Zones::xmax};
			if (z.get_face(zone, face) < 0) continue;
			for (int j {}; j < ny; ++j)
			{
				EXPECT_DOUBLE_EQ(z.z_shift(zone, i, j),
					g.z_shift(origin[zone][0] + i, origin[zone][1] + j));
			}
		}
	}
	EXPECT_DOUBLE_EQ(z.z_period(0), period);
	EXPECT_NEAR(z.wrap_z(0, 3.5 * period + 0.25), 0.5 * period + 0.25,
		1e-12);
	EXPECT_NEAR(z.wrap_z(0, -2.0 * period + 0.1), 0.1, 1e-12);
}

// The twist zShift gives per y cell is the field line's,
// (B . grad z) / (B . grad y), from the background's own field and
// reciprocal basis
TEST_F(FieldAlignedTest, TwistMatchesBasis)
{
	const auto bkg {read()};
	const Zones::Zones& z {bkg.get_zones()};
	for (int zone {}; zone < 6; ++zone)
	for (int i {}; i < z.get_nx(zone); ++i)
	for (int j {}; j < z.get_ny(zone); ++j)
	for (int k : {0, 7, nz - 1})
	{
		// Past a target the ghost cells copy the edge, so there's no
		// twist to compare with
		if ((j == 0 && z.get_face(zone, Zones::ymin) < 0)
			|| (j == z.get_ny(zone) - 1 && z.get_face(zone, Zones::ymax) < 0))
		{
			continue;
		}
		const std::array<double, 3> B {bkg.get_bX()(zone, 0, i, j, k),
			bkg.get_bY()(zone, 0, i, j, k), bkg.get_bZ()(zone, 0, i, j, k)};
		const double B_grad_y {B[0] * bkg.get_dydX()(zone, i, j, k)
			+ B[1] * bkg.get_dydY()(zone, i, j, k)
			+ B[2] * bkg.get_dydZ()(zone, i, j, k)};
		const double B_grad_z {B[0] * bkg.get_dzdX()(zone, i, j, k)
			+ B[1] * bkg.get_dzdY()(zone, i, j, k)
			+ B[2] * bkg.get_dzdZ()(zone, i, j, k)};
		const double twist {0.5 * (z.z_shift(zone, i, j + 1)
			- z.z_shift(zone, i, j - 1))};
		EXPECT_NEAR(B_grad_z / B_grad_y, twist, 1e-10 * std::abs(twist));
	}
}

// The spectral shift is exact for resolved modes, a whole-cell shift is a
// roll, and the Nyquist mode of an even row stays real
TEST(ShiftPeriodic, Exact)
{
	for (int n : {7, 8, 32, 81})
	{
		const double two_pi {2.0 * std::numbers::pi};
		auto f = [&](double k) {return 0.3 + std::sin(two_pi * k / n)
			- 0.7 * std::cos(two_pi * 3 * (k - 0.2) / n);};
		std::vector<double> row(n);
		for (int k {}; k < n; ++k) row[k] = f(k);
		for (double shift : {0.37, -5.9, 2.0 * n + 1.25})
		{
			const auto out {Axisym::shift_periodic(row, shift)};
			for (int k {}; k < n; ++k)
				EXPECT_NEAR(out[k], f(k + shift), 1e-12) << n << " " << shift;
		}
		const auto rolled {Axisym::shift_periodic(row, 2.0)};
		for (int k {}; k < n; ++k)
			EXPECT_NEAR(rolled[k], row[(k + 2) % n], 1e-12);
	}

	// Nyquist mode of an even row: cos(pi k), shifted by half a cell, is 0
	std::vector<double> nyq(8);
	for (int k {}; k < 8; ++k) nyq[k] = k % 2 ? -1.0 : 1.0;
	for (double v : Axisym::shift_periodic(nyq, 0.5))
		EXPECT_NEAR(v, 0.0, 1e-12);
}

// Plasma constant along the field line is interpolated along it, with
// only the error of linear interpolation in z, including across the
// core's branch cut, from one PFR leg to the other and across the wedge's
// toroidal ends. At fixed z the same points are badly off.
TEST_F(FieldAlignedTest, InterpolateScalar)
{
	const auto bkg {read()};
	const Zones::Zones& z {bkg.get_zones()};
	const double ne0 {1e19};

	// The field line's value at a point: in each x column, zShift at the
	// point's y, interpolated between the bracketing planes
	auto expected = [&](int zone, double x, double y, double zp)
	{
		const auto [ix0, ix1] = z.bracket(zone, 0, x);
		const auto [iy0, iy1] = z.bracket(zone, 1, y);
		const double xc0 {z.center(zone, 0, ix0)};
		const double xc1 {z.center(zone, 0, ix1)};
		const double yc0 {z.center(zone, 1, iy0)};
		const double yc1 {z.center(zone, 1, iy1)};
		const double wx {std::clamp((x - xc0) / (xc1 - xc0), 0.0, 1.0)};
		const double wy {std::clamp((y - yc0) / (yc1 - yc0), 0.0, 1.0)};
		double out {};
		for (int a {}; a < 2; ++a)
		{
			const int ix {a ? ix1 : ix0};
			const double zs {(1 - wy) * z.z_shift(zone, ix, iy0)
				+ wy * z.z_shift(zone, ix, iy1)};
			out += (a ? wx : 1 - wx) * (1.0 + 0.5 * profile(zp - zs));
		}
		return out * ne0;
	};

	// Points: (zone, x, y, z) in each zone's coordinates (global grid
	// indices; cell centers are at index + 0.5). None is next to a corner
	// ghost cell at the X-point, whose value is the average of cells on
	// different field lines.
	const std::vector<std::array<double, 4>> points {
		{1, 3.3, 5.8, 0.40},	// core, across the branch cut at the top
		{1, 3.7, 2.2, 1.50},	// core, across the cut at the bottom
		{1, 2.9, 4.1, 0.01},	// core, near the start of the wedge
		{0, 4.1, 1.9, 0.80},	// inner PFR, toward the outer PFR
		{2, 2.6, 6.1, 1.20},	// outer PFR, toward the inner PFR
		{4, 6.4, 3.6, 1.55},	// main SOL, near the end of the wedge
		{3, 7.2, 1.95, 0.70},	// inner leg, toward the SOL
		{5, 6.3, 6.4, 0.30},	// outer leg, from the SOL
		{4, 5.1, 4.5, 0.90}};	// SOL, at the separatrix (x ghost)
	double worst_fixed {};
	for (const auto& [zone_d, x, y, zp] : points)
	{
		const int zone {static_cast<int>(zone_d)};
		for (double t : {0.0, 3e-9})
		{
			const double want {expected(zone, x, y, zp)};
			const double got {bkg.interp_ne(zone, t, x, y, zp)};

			// Linear interpolation in z of a cosine with 32 points per
			// period: at most (2 pi / 32)^2 / 8 of the amplitude
			EXPECT_NEAR(got, want, 0.5 * ne0 * 6e-3) << zone << " " << x
				<< " " << y << " " << zp;
			const double fixed {bkg.interp_4d(bkg.get_ne(), zone, t, x, y,
				zp)};
			worst_fixed = std::max(worst_fixed, std::abs(fixed - want));
		}
	}
	EXPECT_GT(worst_fixed, 0.5 * ne0 * 0.2);
}

// A field that doesn't vary in z interpolates the same either way
TEST_F(FieldAlignedTest, UniformInZ)
{
	const auto bkg {read()};
	const Zones::Zones& z {bkg.get_zones()};
	auto f {Axisym::zoned(z, 2)};
	for (int zone {}; zone < 6; ++zone)
	for (int t {}; t < 2; ++t)
	for (int i {}; i < z.get_nx(zone); ++i)
	for (int j {}; j < z.get_ny(zone); ++j)
	for (int k {}; k < nz; ++k)
	{
		f(zone, t, i, j, k) = 1.0 + 0.3 * (origin[zone][0] + i)
			+ 0.07 * (origin[zone][1] + j) * (origin[zone][1] + j) + t;
	}
	z.fill_ghosts(f);
	for (int zone {}; zone < 6; ++zone)
	for (double x : {0.1, 0.5, 0.97})
	for (double y : {0.0, 0.3, 0.99})
	for (double zp : {0.0, 0.77, period - 1e-3})
	{
		const double xp {z.get_grid_x(zone).front() + x
			* (z.get_grid_x(zone).back() - z.get_grid_x(zone).front())};
		const double yp {z.get_grid_y(zone).front() + y
			* (z.get_grid_y(zone).back() - z.get_grid_y(zone).front())};
		const auto& times {bkg.get_times()};
		const double t {0.4 * times.front() + 0.6 * times.back()};
		EXPECT_NEAR(Utilities::interp_vec4d(f, zone, z, times, t, xp, yp,
			zp, true), Utilities::interp_vec4d(f, zone, z, times, t, xp, yp,
			zp, false), 1e-12);
	}
}

// A vector field that is the same in the (R, phi, Z) frame everywhere
// comes back exactly, whatever the shift: each sample is rotated from
// its cell's toroidal angle to the point's, including samples from the
// periodic ghost cells at the ends of the wedge (whose vectors are at the
// ghost's own angle)
TEST_F(FieldAlignedTest, InterpolateVector)
{
	const auto bkg {read()};
	const Zones::Zones& z {bkg.get_zones()};
	auto vX {Axisym::zoned(z, 2)};
	auto vY {Axisym::zoned(z, 2)};
	auto vZ {Axisym::zoned(z, 2)};
	const double vR {0.3}, vphi {-2.0}, vz {0.6};
	for (int zone {}; zone < 6; ++zone)
	for (int t {}; t < 2; ++t)
	for (int i {}; i < z.get_nx(zone); ++i)
	for (int j {}; j < z.get_ny(zone); ++j)
	for (int k {}; k < nz; ++k)
	{
		const double phi {z.get_z(zone)[k]};
		vX(zone, t, i, j, k) = vR * std::cos(phi) - vphi * std::sin(phi);
		vY(zone, t, i, j, k) = vR * std::sin(phi) + vphi * std::cos(phi);
		vZ(zone, t, i, j, k) = vz;
	}
	// Ghost cells past the wedge's z faces hold the vector rotated to
	// their own angle, as the background's vectors do
	z.fill_ghosts(vX, vY, bkg.get_phi_sign());
	z.fill_ghosts(vZ);

	for (int zone {}; zone < 6; ++zone)
	for (double y : {0.05, 0.5, 0.95})
	for (double zp : {0.0, 0.02, 0.6, period - 0.01})
	{
		const double xp {z.get_x(zone).front() + 0.3};
		const double yp {z.get_grid_y(zone).front() + y
			* (z.get_grid_y(zone).back() - z.get_grid_y(zone).front())};
		const auto v {bkg.interp_vec3(vX, vY, vZ, zone, bkg.get_times()[0],
			xp, yp, zp)};
		EXPECT_NEAR(v[0], vR * std::cos(zp) - vphi * std::sin(zp), 1e-12);
		EXPECT_NEAR(v[1], vR * std::sin(zp) + vphi * std::cos(zp), 1e-12);
		EXPECT_NEAR(v[2], vz, 1e-12);
	}
}

// The y derivative of a field-aligned quantity: along the field line the
// shifted neighboring planes match the cell's exactly, so df/dy|_zeta is
// 0, and at fixed z that leaves -(dzShift/dy) df/dz, which is the
// analytic derivative up to the central difference in z
TEST_F(FieldAlignedTest, AlignedDerivative)
{
	const auto bkg {read()};
	const Zones::Zones& z {bkg.get_zones()};
	const auto& vp {bkg.get_vp()};
	const auto dfdy {Axisym::aligned_dfdy(z, vp)};
	const double dz {period / nz};
	const double two_pi {2.0 * std::numbers::pi};
	for (int zone {}; zone < 6; ++zone)
	for (int i {}; i < z.get_nx(zone); ++i)
	for (int j {}; j < z.get_ny(zone); ++j)
	for (int k {}; k < nz; ++k)
	{
		EXPECT_NEAR(dfdy(zone, 1, i, j, k), 0.0, 1e-9 * tnorm);

		const int gi {origin[zone][0] + i};
		const double s {m_grid.twist_per_cell(gi)};
		const auto d {Axisym::derivatives(z, vp, zone, 1, i, j, k, &dfdy)};
		EXPECT_NEAR(d[1], -s * d[2], 1e-9 * tnorm);

		const double zeta {k * dz - z.z_shift(zone, i, j)};
		const double analytic {0.3 * tnorm * s * two_pi / period
			* std::sin(two_pi * zeta / period)};
		EXPECT_NEAR(d[1], analytic, 0.02 * tnorm * std::abs(s) * two_pi
			/ period);
	}
}

// The parallel gradient of a quantity that grows along the field line,
// f = y p(zeta): df/dy|_zeta = p(zeta) exactly, and
// b . grad f = (B . grad y / B) p(zeta) = (Bz / (dZ B)) p(zeta) on the
// test grid. Checked away from where y jumps between zones (the core's
// branch cut and the PFR), where y p isn't continuous along the field.
TEST_F(FieldAlignedTest, GradPar)
{
	const auto bkg {read()};
	const Zones::Zones& z {bkg.get_zones()};
	const double dz {period / nz};
	auto f {Axisym::zoned(z, 2)};
	for (int zone {}; zone < 6; ++zone)
	for (int t {}; t < 2; ++t)
	for (int i {}; i < z.get_nx(zone); ++i)
	for (int j {}; j < z.get_ny(zone); ++j)
	for (int k {}; k < nz; ++k)
	{
		f(zone, t, i, j, k) = z.get_y(zone)[j] * profile(k * dz
			- z.z_shift(zone, i, j));
	}
	z.fill_ghosts(f);
	const auto dfdy {Axisym::aligned_dfdy(z, f)};
	for (int zone : {1, 4})
	for (int i {}; i < z.get_nx(zone); ++i)
	for (int j {1}; j < z.get_ny(zone) - 1; ++j)
	for (int k {}; k < nz; ++k)
	{
		const double p {profile(k * dz - z.z_shift(zone, i, j))};
		const double Bt {m_grid.Bt_at(origin[zone][0] + i)};
		const double B {std::hypot(m_grid.Bz, Bt)};
		EXPECT_NEAR(Axisym::grad_par(bkg, dfdy, zone, 0, i, j, k),
			m_grid.Bz / (m_grid.dZ * B) * p, 1e-10);
	}
}

// E from a potential that is constant along the field line has no
// parallel component, in every cell including next to the core's branch
// cut, the PFR's connection and the targets
TEST_F(FieldAlignedTest, NoParallelE)
{
	const auto bkg {read()};
	const Zones::Zones& z {bkg.get_zones()};
	double max_E {}, max_par {};
	for (int zone {}; zone < 6; ++zone)
	for (int t {}; t < 2; ++t)
	for (int i {}; i < z.get_nx(zone); ++i)
	for (int j {}; j < z.get_ny(zone); ++j)
	for (int k {}; k < nz; ++k)
	{
		const std::array<double, 3> E {bkg.get_eX()(zone, t, i, j, k),
			bkg.get_eY()(zone, t, i, j, k), bkg.get_eZ()(zone, t, i, j, k)};
		const std::array<double, 3> B {bkg.get_bX()(zone, t, i, j, k),
			bkg.get_bY()(zone, t, i, j, k), bkg.get_bZ()(zone, t, i, j, k)};
		const double Bmag {std::sqrt(B[0] * B[0] + B[1] * B[1]
			+ B[2] * B[2])};
		max_E = std::max(max_E, std::sqrt(E[0] * E[0] + E[1] * E[1]
			+ E[2] * E[2]));
		max_par = std::max(max_par, std::abs(E[0] * B[0] + E[1] * B[1]
			+ E[2] * B[2]) / Bmag);
	}

	// E isn't trivially zero (about 2 pi / period phi |grad zeta|), but
	// its parallel part is, to rounding
	EXPECT_GT(max_E, 10.0);
	EXPECT_LT(max_par, 1e-12 * max_E);
}

namespace
{
	// Electron and ion q_par expected from kappa (Hermes-3's
	// "multispecies" collision times, ni = ne) and a parallel gradient
	std::array<double, 2> expected_q(double ne, double te, double ti,
		double gpar_e, double gpar_i)
	{
		const double mi {2.0 * Braginskii::mp};
		const double nu_ei {Braginskii::nu_ei(ne, ne, te, ti, 1.0, 2.0)};
		const double nu_e {Braginskii::nu_ee(ne, te) + nu_ei};
		const double nu_i {Braginskii::nu_ii(ne, ti, 1.0, 2.0, ne, ti, 1.0,
			2.0) + nu_ei * Braginskii::me / mi};
		return {-Braginskii::kappa_par(Braginskii::kappa_coef_e_default, ne,
			te, nu_e, Braginskii::me) * gpar_e,
			-Braginskii::kappa_par(Braginskii::kappa_coef_i_default, ne, ti,
			nu_i, mi) * gpar_i};
	}
}

// A HERMES background gets both heat fluxes by default, and none with
// hermes_heat_flux = "off"
TEST_F(FieldAlignedTest, HeatFluxOption)
{
	const auto bkg {read()};
	EXPECT_FALSE(bkg.get_q_par_e().get_data().empty());
	EXPECT_FALSE(bkg.get_q_par_i().get_data().empty());

	Options::Options opts {};
	opts.set_bkg_source("hermes");
	opts.set_hermes_grid_file((m_dir / "grid.nc").string());
	opts.set_hermes_snap_dir((m_dir / "snaps").string());
	opts.set_hermes_snap_end(-1);
	opts.set_hermes_tnorm(tnorm);
	opts.set_hermes_zperiod(zperiod);
	opts.set_hermes_heat_flux("off");
	const auto off {Hermes::read_hermes(opts)};
	EXPECT_TRUE(off.get_q_par_e().get_data().empty());
	EXPECT_TRUE(off.get_q_par_i().get_data().empty());
}

// With T rising along the field line, q_par = -kappa grad_par T for both
// species, so it points back down the gradient. Here T = T0 + y (1 +
// p(zeta) / 2) for the electrons and twice that rise for the ions, whose
// parallel gradient is (B . grad y / B) (1 + p / 2) (see GradPar).
TEST_F(FieldAlignedTest, HeatFluxAlongField)
{
	auto bkg {read()};
	const Zones::Zones& z {bkg.get_zones()};
	const double dz {period / nz};
	auto& te {bkg.get_te()};
	auto& ti {bkg.get_ti()};
	for (int zone {}; zone < 6; ++zone)
	for (int t {}; t < 2; ++t)
	for (int i {}; i < z.get_nx(zone); ++i)
	for (int j {}; j < z.get_ny(zone); ++j)
	for (int k {}; k < nz; ++k)
	{
		const double rise {z.get_y(zone)[j] * (1.0 + 0.5 * profile(k * dz
			- z.z_shift(zone, i, j)))};
		te(zone, t, i, j, k) = tnorm + rise;
		ti(zone, t, i, j, k) = tnorm + 2.0 * rise;
	}
	z.fill_ghosts(te);
	z.fill_ghosts(ti);
	Axisym::calc_conductive_heat_fluxes(bkg, Axisym::Conduction {});

	const auto& ne {bkg.get_ne()};
	for (int zone : {1, 4})
	for (int i {}; i < z.get_nx(zone); ++i)
	for (int j {1}; j < z.get_ny(zone) - 1; ++j)
	for (int k {}; k < nz; ++k)
	{
		const double p {profile(k * dz - z.z_shift(zone, i, j))};
		const double Bt {m_grid.Bt_at(origin[zone][0] + i)};
		const double G {m_grid.Bz / (m_grid.dZ * std::hypot(m_grid.Bz, Bt))
			* (1.0 + 0.5 * p)};
		const auto q {expected_q(ne(zone, 0, i, j, k), te(zone, 0, i, j, k),
			ti(zone, 0, i, j, k), G, 2.0 * G)};
		EXPECT_NEAR(bkg.get_q_par_e()(zone, 0, i, j, k), q[0],
			1e-9 * std::abs(q[0]));
		EXPECT_NEAR(bkg.get_q_par_i()(zone, 0, i, j, k), q[1],
			1e-9 * std::abs(q[1]));

		// Down the gradient: B . grad y > 0 here, so T rises along B
		EXPECT_LT(bkg.get_q_par_e()(zone, 0, i, j, k), 0.0);
	}
}

// T that varies only across the field (constant along field lines) has
// no parallel heat flux, although its gradient is large. b . grad T
// would pick up part of it (see grad_par).
TEST_F(FieldAlignedTest, NoHeatFluxAcrossField)
{
	auto bkg {read()};
	const Zones::Zones& z {bkg.get_zones()};
	const double dz {period / nz};
	auto& te {bkg.get_te()};
	auto& ti {bkg.get_ti()};
	for (int zone {}; zone < 6; ++zone)
	for (int t {}; t < 2; ++t)
	for (int i {}; i < z.get_nx(zone); ++i)
	for (int j {}; j < z.get_ny(zone); ++j)
	for (int k {}; k < nz; ++k)
	{
		const double p {profile(k * dz - z.z_shift(zone, i, j))};
		te(zone, t, i, j, k) = tnorm * (1.0 + 0.3 * p);
		ti(zone, t, i, j, k) = tnorm * (1.0 + 0.3 * p);
	}
	z.fill_ghosts(te);
	z.fill_ghosts(ti);
	Axisym::calc_conductive_heat_fluxes(bkg, Axisym::Conduction {});

	// Scale: the heat flux if this T's toroidal gradient were along B
	const double kappa_e {-expected_q(1e19, tnorm, tnorm, 1.0, 1.0)[0]};
	const double grad {0.3 * tnorm * 2.0 * std::numbers::pi / period
		/ m_grid.R0};
	double max_q {};
	for (const double q : bkg.get_q_par_e().get_data())
		max_q = std::max(max_q, std::abs(q));
	EXPECT_LT(max_q, 1e-9 * kappa_e * grad);
}
