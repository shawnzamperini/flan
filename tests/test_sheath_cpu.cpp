/**
* @file test_sheath_cpu.cpp
*
* @brief Tests of the sheath field (sheath.cpp): the split of the drop,
*        both potential profiles, and the field at particles in front of a
*        slab wall
*/
#include <array>
#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "background.h"
#include "boris.h"
#include "indices.h"
#include "options.h"
#include "read_bkg.h"
#include "slots.h"
#include "sheath.h"
#include "sputtering.h"

namespace
{
	constexpr double pi {3.141592653589793};
	constexpr double me_amu {5.4858e-4};
	constexpr double mi_amu {2.014};

	// A D plasma at Te = Ti = 20 eV, 1e19 m^-3 and 2 T, with the field at
	// deg degrees from the wall's normal
	Sheath::Params params_at(const double deg, const double ti = 20.0)
	{
		return Sheath::make_params(20.0, ti, 1e19, 2.0,
			std::cos(deg * pi / 180.0), me_amu, mi_amu);
	}

	// The integral of E . n from the wall out to the field's range, by
	// the midpoint rule on a grid fine on the shortest scale (lambda_D)
	double integrate_field(const Sheath::Params& p, const Sheath::Model m)
	{
		const double r {Sheath::range(p, m)};
		const double h {0.01 * p.lambda_d};
		const int n {static_cast<int>(std::ceil(r / h))};
		const double dd {r / n};
		double sum {};
		for (int i {}; i < n; ++i)
			sum += Sheath::field(p, m, (i + 0.5) * dd) * dd;
		return sum;
	}

	const std::vector<double> angles {0.0, 10.0, 30.0, 45.0, 60.0, 75.0,
		80.0, 85.0, 86.0, 87.0, 89.0};
}

// The floating drop, the Debye length and the sound gyroradius, by hand
TEST(SheathParams, Scales)
{
	const auto p {params_at(45.0)};
	const double lam_w {0.5 * std::log(2.0 * pi * me_amu / mi_amu * 2.0)};
	EXPECT_NEAR(p.lam_w, lam_w, 1e-12);
	EXPECT_NEAR(p.lam_w, -2.839, 1e-3);
	EXPECT_NEAR(p.phi_wall, 20.0 * lam_w, 1e-10);
	EXPECT_NEAR(p.lambda_d, std::sqrt(8.85e-12 * 20.0 / (1e19 * 1.602e-19)),
		1e-15);
	const double mi {mi_amu * 1.661e-27};
	const double c_s {std::sqrt(1.602e-19 * 40.0 / mi)};
	EXPECT_NEAR(p.rho_s, c_s / (1.602e-19 * 2.0 / mi), 1e-12);
	EXPECT_NEAR(p.theta, pi / 4.0, 1e-12);
}

// f_D = 1 - ln cos(theta) / lam_w, 1 at normal incidence, falling to 0 at
// the critical angle (about 86.6 degrees for D at Ti = Te) and 0 past it
TEST(SheathParams, DebyeShare)
{
	for (const double deg : angles)
	{
		const auto p {params_at(deg)};
		const double expect {std::max(1.0 - std::log(std::cos(deg * pi
			/ 180.0)) / p.lam_w, 0.0)};
		EXPECT_NEAR(p.f_d, expect, 1e-12) << deg;
	}
	EXPECT_DOUBLE_EQ(params_at(0.0).f_d, 1.0);
	const double crit {std::acos(std::exp(params_at(0.0).lam_w)) * 180.0
		/ pi};
	EXPECT_NEAR(crit, 86.65, 0.01);
	EXPECT_GT(params_at(crit - 0.1).f_d, 0.0);
	EXPECT_DOUBLE_EQ(params_at(crit + 0.1).f_d, 0.0);
	EXPECT_FALSE(params_at(crit + 0.1).debye);

	// Without a field there's no presheath
	const auto p {Sheath::make_params(20.0, 20.0, 1e19, 0.0, 0.0, me_amu,
		mi_amu)};
	EXPECT_DOUBLE_EQ(p.f_d, 1.0);

	// No plasma, no sheath
	const auto none {Sheath::make_params(0.0, 20.0, 1e19, 2.0, 1.0, me_amu,
		mi_amu)};
	EXPECT_DOUBLE_EQ(none.phi_wall, 0.0);
	EXPECT_DOUBLE_EQ(Sheath::field(none, Sheath::Model::brooks, 0.0), 0.0);
}

// Both profiles are phi_wall at the wall and 0 far away, and the field
// integrates to the drop: int_0^inf E . n dd = phi(0) - phi(inf)
TEST(SheathProfile, IntegratesToDrop)
{
	for (const auto m : {Sheath::Model::brooks, Sheath::Model::borodkina})
	for (const double deg : angles)
	{
		const auto p {params_at(deg)};
		EXPECT_NEAR(Sheath::potential(p, m, 0.0), p.phi_wall, 1e-10)
			<< deg;
		EXPECT_NEAR(Sheath::potential(p, m, Sheath::range(p, m)), 0.0,
			1e-9 * std::abs(p.phi_wall)) << deg;
		EXPECT_NEAR(integrate_field(p, m), p.phi_wall,
			1e-4 * std::abs(p.phi_wall)) << static_cast<int>(m) << " "
			<< deg;
	}
}

// The field is minus the potential's slope, and points into the wall
// (E . n < 0) everywhere
TEST(SheathProfile, FieldIsSlope)
{
	for (const auto m : {Sheath::Model::brooks, Sheath::Model::borodkina})
	for (const double deg : angles)
	{
		const auto p {params_at(deg)};
		const double r {Sheath::range(p, m, 7.0)};
		for (int i {}; i <= 50; ++i)
		{
			const double d {r * i / 50.0};
			const double h {1e-4 * p.lambda_d};
			const double slope {(Sheath::potential(p, m, d + h)
				- Sheath::potential(p, m, std::max(d - h, 0.0)))
				/ (d + h - std::max(d - h, 0.0))};
			const double e {Sheath::field(p, m, d)};
			EXPECT_LT(e, 0.0) << deg << " " << d;
			EXPECT_NEAR(e, -slope, 1e-4 * std::abs(e) + 1e-6) << deg << " "
				<< d;
		}
	}
}

// Brooks' profile: the Debye part decays over 2 lambda_D and holds f_D of
// the drop, the presheath part over rho_s
TEST(SheathProfile, Brooks)
{
	const auto p {params_at(60.0)};
	const auto m {Sheath::Model::brooks};
	const double d {3.0 * p.lambda_d};
	const double expect {p.phi_wall * (p.f_d * std::exp(-1.5)
		+ (1.0 - p.f_d) * std::exp(-d / p.rho_s))};
	EXPECT_NEAR(Sheath::potential(p, m, d), expect, 1e-12);

	// Far past the Debye sheath, only the presheath part is left
	const double far {50.0 * p.lambda_d};
	EXPECT_NEAR(Sheath::potential(p, m, far), p.phi_wall * (1.0 - p.f_d)
		* std::exp(-far / p.rho_s), 1e-9);
}

// Borodkina and Tsvetkov's profile: lam_mps = ln cos(theta) at xi_mps,
// with the potential and field continuous there, the presheath decaying
// over rho_s sin(theta), and a single Debye layer at normal incidence and
// a single presheath past the critical angle
TEST(SheathProfile, Borodkina)
{
	const auto m {Sheath::Model::borodkina};
	for (const double deg : {10.0, 30.0, 45.0, 60.0, 75.0, 80.0, 85.0})
	{
		const auto p {params_at(deg)};
		ASSERT_TRUE(p.debye) << deg;
		EXPECT_GT(p.xi_mps, 0.0);
		EXPECT_NEAR(p.lam_mps, std::log(std::cos(deg * pi / 180.0)),
			1e-12);
		EXPECT_NEAR(p.l_mps * p.lambda_d, p.rho_s * std::sin(p.theta),
			1e-12);
		const double d0 {p.xi_mps * p.lambda_d};
		const double in {d0 * (1.0 - 1e-10)};
		const double out {d0 * (1.0 + 1e-10)};
		EXPECT_NEAR(Sheath::potential(p, m, in), p.te * p.lam_mps, 1e-6)
			<< deg;
		EXPECT_NEAR(Sheath::potential(p, m, out), p.te * p.lam_mps, 1e-6)
			<< deg;
		const double e_in {Sheath::field(p, m, in)};
		EXPECT_NEAR(Sheath::field(p, m, out), e_in, 1e-3 * std::abs(e_in))
			<< deg;
	}

	// At normal incidence: lam_w + Q (1 - exp(-a xi)) all the way, with
	// Q = -lam_w and 1 / a close to Brooks' 2 lambda_D
	const auto p0 {params_at(0.0)};
	ASSERT_TRUE(p0.debye);
	EXPECT_FALSE(std::isfinite(p0.xi_mps));
	EXPECT_NEAR(p0.q, -p0.lam_w, 1e-12);
	EXPECT_NEAR(1.0 / p0.a, 2.0, 0.1);

	// Past the critical angle: lam_w exp(-xi / L)
	const auto p89 {params_at(89.0)};
	ASSERT_FALSE(p89.debye);
	const double d {p89.rho_s};
	EXPECT_NEAR(Sheath::potential(p89, m, d), p89.phi_wall
		* std::exp(-1.0 / std::sin(p89.theta)), 1e-9);
}

// sheath_fd and sheath_drop: f_D and the drop given, as in Brooks'
// reference case (f_D = 0.25, 3 Te at 87 degrees, past the critical
// angle). Both profiles still integrate to the drop, and Borodkina's is
// still continuous at xi_mps.
TEST(SheathProfile, Overrides)
{
	const auto p {Sheath::make_params(75.0, 75.0, 1e20, 5.0,
		std::cos(87.0 * pi / 180.0), me_amu, 2.5, 0.25, 3.0)};
	EXPECT_DOUBLE_EQ(p.f_d, 0.25);
	EXPECT_DOUBLE_EQ(p.lam_w, -3.0);
	EXPECT_DOUBLE_EQ(p.phi_wall, -225.0);
	EXPECT_DOUBLE_EQ(p.lam_mps, -2.25);
	ASSERT_TRUE(p.debye);
	for (const auto m : {Sheath::Model::brooks, Sheath::Model::borodkina})
	{
		EXPECT_NEAR(Sheath::potential(p, m, 0.0), -225.0, 1e-9);
		EXPECT_NEAR(integrate_field(p, m), -225.0, 1e-4 * 225.0);
	}
	const double d0 {p.xi_mps * p.lambda_d};
	const auto m {Sheath::Model::borodkina};
	EXPECT_NEAR(Sheath::potential(p, m, d0 * (1.0 - 1e-10)),
		Sheath::potential(p, m, d0 * (1.0 + 1e-10)), 1e-6);

	// Brooks: the Debye part holds f_D of the drop
	const auto b {Sheath::Model::brooks};
	EXPECT_NEAR(Sheath::potential(p, b, 10.0 * p.lambda_d), -225.0 * (0.25
		* std::exp(-5.0) + 0.75 * std::exp(-10.0 * p.lambda_d / p.rho_s)),
		1e-9);

	// Only the drop given: f_D from the angle, of the new drop
	const auto q {Sheath::make_params(20.0, 20.0, 1e19, 2.0,
		std::cos(60.0 * pi / 180.0), me_amu, mi_amu, -1.0, 4.0)};
	EXPECT_NEAR(q.f_d, 1.0 - std::log(0.5) / -4.0, 1e-12);
}

// "off" gives nothing
TEST(SheathProfile, Off)
{
	const auto p {params_at(45.0)};
	EXPECT_DOUBLE_EQ(Sheath::potential(p, Sheath::Model::off, 0.0), 0.0);
	EXPECT_DOUBLE_EQ(Sheath::field(p, Sheath::Model::off, 0.0), 0.0);
	EXPECT_DOUBLE_EQ(Sheath::range(p, Sheath::Model::off), 0.0);
}

namespace
{
	// The 5 cm slab with a wall at y min, 8 cells a side, and the field
	// at b_angle degrees to the wall
	Background::Background slab_wall(const double b_angle)
	{
		Options::Options opts {};
		opts.set_imp_atom_num(74);
		opts.set_bkg_source("slab");
		opts.set_slab_ne(1e19);
		opts.set_slab_te(20.0);
		opts.set_slab_ti(20.0);
		opts.set_slab_bmag(2.0);
		opts.set_slab_tmax(1e-5);
		opts.set_slab_wall("ymin");
		opts.set_slab_b_angle(b_angle);
		return Background::read_bkg(opts);
	}
}

// In front of a slab wall at y min (normal +Y, at Y = -2.5 cm) the field
// is E . n at the particle's height above it, along -Y, in the cells next
// to the wall only. With the field at 30 degrees to the wall, theta is 60
// degrees.
TEST(SheathField, SlabWall)
{
	const auto bkg {slab_wall(30.0)};
	const Sputtering::Surface surface(bkg);
	for (const auto m : {Sheath::Model::brooks, Sheath::Model::borodkina})
	{
		const Sheath::Sheath sheath(surface, bkg, m);
		ASSERT_EQ(sheath.size(), 64u);
		EXPECT_EQ(sheath.nframes(), static_cast<int>(bkg.get_times().size()));
		EXPECT_EQ(sheath.num_thin(), 0);

		const auto& p {sheath.params(0, 0)};
		EXPECT_NEAR(p.theta * 180.0 / pi, 60.0, 1e-9);
		const double lam_w {0.5 * std::log(2.0 * pi
			* bkg.get_elec_mass_amu() / bkg.get_ion_mass_amu() * 2.0)};
		EXPECT_NEAR(p.lam_w, lam_w, 1e-12);
		EXPECT_NEAR(p.f_d, 1.0 - std::log(0.5) / lam_w, 1e-12);

		const double side {0.05 / 8.0};
		for (const double d : {0.0, 1e-6, 1e-5, 1e-4, 1e-3, 5e-3})
		{
			// x in cell 3, z in cell 5: y is d above the wall
			const std::array<double, 3> pos {3.5 * side, -0.025 + d,
				-0.025 + 5.5 * side};
			std::size_t e {};
			double dist {};
			ASSERT_TRUE(sheath.distance(0, 3, 0, 5, pos, e, dist));
			EXPECT_NEAR(dist, d, 1e-15);
			EXPECT_EQ(surface.elements()[e].i, 3);
			EXPECT_EQ(surface.elements()[e].k, 5);

			const auto E {sheath.field_at(0, 0, 3, 0, 5, pos)};
			const double en {d <= Sheath::range(p, m)
				? Sheath::field(p, m, d) : 0.0};
			EXPECT_DOUBLE_EQ(E[0], 0.0);
			EXPECT_NEAR(E[1], en, 1e-9 * std::abs(en) + 1e-12) << d;
			EXPECT_DOUBLE_EQ(E[2], 0.0);
			if (d < 1e-3) EXPECT_LT(E[1], 0.0) << d;
		}

		// The next row of cells isn't next to the wall: no field
		const std::array<double, 3> pos {3.5 * side, -0.025 + 1.5 * side,
			0.0};
		std::size_t e {};
		double dist {};
		EXPECT_FALSE(sheath.distance(0, 3, 1, 5, pos, e, dist));
		const auto E {sheath.field_at(0, 0, 3, 1, 5, pos)};
		for (const double c : E) EXPECT_DOUBLE_EQ(c, 0.0);
	}

	// Off, nothing at all
	const Sheath::Sheath off(surface, bkg, Sheath::Model::off);
	const auto E {off.field_at(0, 0, 3, 0, 5, {0.02, -0.025, 0.0})};
	for (const double c : E) EXPECT_DOUBLE_EQ(c, 0.0);
}

// A positive ion released at rest in front of the wall is pulled into it,
// gaining Z e (phi(start) - phi_wall) on the way: integrate its motion along the
// normal through the field alone (no B), with small steps
TEST(SheathField, IonGainsDrop)
{
	const auto bkg {slab_wall(30.0)};
	const Sputtering::Surface surface(bkg);
	for (const auto m : {Sheath::Model::brooks, Sheath::Model::borodkina})
	{
		const Sheath::Sheath sheath(surface, bkg, m);
		const auto& p {sheath.params(0, 0)};
		const double start {Sheath::range(p, m, 7.0)};
		ASSERT_LT(start, 0.05 / 8.0);
		const double gain {Sheath::potential(p, m, start) - p.phi_wall};

		// W+ (the mass doesn't change the energy)
		const double q_m {1.602e-19 / (183.84 * 1.661e-27)};
		double y {start};
		double v {};
		for (int n {}; n < 10000000 && y > 0.0; ++n)
		{
			// At most 0.002 lambda_D a step, so the Debye layer is
			// resolved
			const double dt {std::min(1e-10, 2e-3 * p.lambda_d
				/ std::max(std::abs(v), 1.0))};
			const auto E {sheath.field_at(0, 0, 3, 0, 5,
				{0.02, -0.025 + y, 0.0})};
			v += q_m * E[1] * dt;
			y += v * dt;
		}
		ASSERT_LE(y, 0.0);
		const double energy_ev {0.5 * v * v / q_m};
		EXPECT_NEAR(energy_ev, gain, 1e-3 * gain)
			<< static_cast<int>(m);
	}
}

namespace
{
	/**
	* Follows one W ion of charge z, released at rest `start` from the slab
	* wall at y min, until it crosses the wall, as the main loop does:
	* each step, Sheath::substeps, the Boris push of the particles that
	* aren't sub-stepped, Boris::substep_cpu, the position step for the
	* rest, and the cell indices found again. Returns its kinetic energy
	* when it crosses (eV), or -1 if it didn't within max_steps.
	*/
	double impact_energy(const Background::Background& bkg,
		const Sheath::Sheath& sheath, const double start, const int z,
		const double dt, const int max_steps, int* max_nsub = nullptr)
	{
		const double mass {183.84 * 1.661e-27};
		Slots::Slots slots(1);
		slots.set_mass(mass);
		slots.set_q(0, z);
		slots.set_x(0, 0.01);
		slots.set_y(0, -0.025 + start);
		slots.set_z(0, 0.0);
		slots.set_vX(0, 0.0);
		slots.set_vY(0, 0.0);
		slots.set_vZ(0, 0.0);
		slots.set_t(0, bkg.get_times()[0]);
		slots.set_tidx(0, 0);
		const auto& gx {bkg.get_grid_x(0)};
		const auto& gy {bkg.get_grid_y(0)};
		const auto& gz {bkg.get_grid_z(0)};
		auto set_cell = [&]()
		{
			slots.set_xidx(0, Indices::get_nearest_cell_index_cpu(gx,
				slots.x()[0]));
			slots.set_yidx(0, Indices::get_nearest_cell_index_cpu(gy,
				slots.y()[0]));
			slots.set_zidx(0, Indices::get_nearest_cell_index_cpu(gz,
				slots.z()[0]));
		};
		set_cell();

		const Options::Options opts {};
		const double q_m {z * 1.602e-19 / mass};
		Boris::update_velocity_cpu(slots, bkg, opts, -0.5 * dt, &sheath);
		std::vector<int> nsub(1, 1);
		for (int n {}; n < max_steps; ++n)
		{
			const double speed {std::hypot(slots.vX()[0], slots.vY()[0],
				slots.vZ()[0])};
			nsub[0] = sheath.substeps(0, 0, slots.xidx()[0],
				slots.yidx()[0], slots.zidx()[0], {slots.x()[0],
				slots.y()[0], slots.z()[0]}, speed, q_m, dt);
			if (max_nsub != nullptr) *max_nsub = std::max(*max_nsub,
				nsub[0]);
			Boris::update_velocity_cpu(slots, bkg, opts, dt, &sheath, &nsub);
			Boris::substep_cpu(slots, bkg, dt, sheath, nsub);
			if (nsub[0] <= 1)
			{
				slots.set_x(0, slots.x()[0] + slots.vx()[0] * dt);
				slots.set_y(0, slots.y()[0] + slots.vy()[0] * dt);
				slots.set_z(0, slots.z()[0] + slots.vz()[0] * dt);
			}
			slots.set_t(0, slots.t()[0] + dt);
			if (slots.y()[0] < gy.front())
			{
				const double v2 {slots.vX()[0] * slots.vX()[0]
					+ slots.vY()[0] * slots.vY()[0] + slots.vZ()[0]
					* slots.vZ()[0]};
				return 0.5 * mass * v2 / 1.602e-19;
			}
			set_cell();
		}
		return -1.0;
	}
}

// Sub-stepping near the wall (S4): an ion released at rest 4 presheath
// lengths from the wall, with B at 30 degrees to it, no
// collisions and no ionization, hits it with z (phi(start) - phi_wall)
// to 1e-3, whatever the time step, from 1e-10 to 1e-8 s. B does no work,
// so that's energy conservation in a static field. (1e-11 s passes too,
// but takes over a minute in a debug build; it uses the same 2 sub-steps
// a step as 1e-10 s.)
TEST(SheathSubstep, ImpactEnergy)
{
	const auto bkg {slab_wall(30.0)};
	const Sputtering::Surface surface(bkg);
	for (const auto m : {Sheath::Model::brooks, Sheath::Model::borodkina})
	{
		const Sheath::Sheath sheath(surface, bkg, m);
		const auto& p {sheath.params(0, 0)};
		const double start {Sheath::range(p, m, 4.0)};
		const double gain {Sheath::potential(p, m, start) - p.phi_wall};
		for (const int z : {1, 4})
		for (const double dt : {1e-10, 1e-9, 1e-8})
		{
			int max_nsub {};
			const double e {impact_energy(bkg, sheath, start, z, dt,
				static_cast<int>(2e-4 / dt) + 10, &max_nsub)};
			ASSERT_GT(e, 0.0) << static_cast<int>(m) << " " << z << " "
				<< dt;
			EXPECT_NEAR(e, z * gain, 1e-3 * z * gain) << static_cast<int>(m)
				<< " Z = " << z << " dt = " << dt << " up to " << max_nsub
				<< " sub-steps";
		}
	}
}

// Without sub-stepping a large step gets the energy badly wrong: what
// S4 fixes. (The same run with substeps forced to 1.)
TEST(SheathSubstep, NeededAtLargeSteps)
{
	const auto bkg {slab_wall(30.0)};
	const Sputtering::Surface surface(bkg);
	const Sheath::Sheath sheath(surface, bkg, Sheath::Model::brooks);
	const auto& p {sheath.params(0, 0)};
	const auto m {Sheath::Model::brooks};
	const double start {Sheath::range(p, m, 7.0)};
	const double gain {Sheath::potential(p, m, start) - p.phi_wall};

	// A plain loop with the sheath field but no sub-steps
	const double mass {183.84 * 1.661e-27};
	const double dt {1e-8};
	Slots::Slots slots(1);
	slots.set_mass(mass);
	slots.set_q(0, 1);
	slots.set_x(0, 0.01);
	slots.set_y(0, -0.025 + start);
	slots.set_t(0, bkg.get_times()[0]);
	const auto& gx {bkg.get_grid_x(0)};
	const auto& gy {bkg.get_grid_y(0)};
	const auto& gz {bkg.get_grid_z(0)};
	const Options::Options opts {};
	double e {-1.0};
	for (int n {}; n < 100000; ++n)
	{
		slots.set_xidx(0, Indices::get_nearest_cell_index_cpu(gx,
			slots.x()[0]));
		slots.set_yidx(0, Indices::get_nearest_cell_index_cpu(gy,
			slots.y()[0]));
		slots.set_zidx(0, Indices::get_nearest_cell_index_cpu(gz,
			slots.z()[0]));
		Boris::update_velocity_cpu(slots, bkg, opts, dt, &sheath);
		slots.set_x(0, slots.x()[0] + slots.vx()[0] * dt);
		slots.set_y(0, slots.y()[0] + slots.vy()[0] * dt);
		slots.set_z(0, slots.z()[0] + slots.vz()[0] * dt);
		if (slots.y()[0] < gy.front())
		{
			e = 0.5 * mass * (slots.vX()[0] * slots.vX()[0] + slots.vY()[0]
				* slots.vY()[0] + slots.vZ()[0] * slots.vZ()[0]) / 1.602e-19;
			break;
		}
	}
	ASSERT_GT(e, 0.0);
	EXPECT_GT(std::abs(e - gain), 0.05 * gain) << e << " vs " << gain;
}

// The impact sampler (sample_impacts): D+ test ions from 10 rho_s through
// either profile, at field angles from normal to past the critical angle.
// Each hits the wall with its starting energy plus e (phi(start) -
// phi_wall), to 1e-3 of the drop (B does no work and the field is
// static), and nearly all of them hit. Weights are the ions' speeds along
// B toward the wall, so they're positive.
TEST(SheathImpacts, EnergyConserved)
{
	for (const auto m : {Sheath::Model::brooks, Sheath::Model::borodkina})
	for (const double deg : {0.0, 30.0, 60.0, 85.0, 89.0})
	{
		const auto p {params_at(deg)};
		const auto r {Sheath::sample_impacts(p, m, 200, 7)};
		EXPECT_EQ(r.num_ions, 200);
		EXPECT_LE(r.lost, 4) << static_cast<int>(m) << " " << deg;
		EXPECT_EQ(static_cast<int>(r.hits.size()) + r.lost, 200);
		EXPECT_LT(r.max_energy_error, 1e-3) << static_cast<int>(m) << " "
			<< deg;
		EXPECT_NEAR(r.start_distance, 10.0 * p.rho_s, 1e-15);
		for (const auto& h : r.hits)
		{
			EXPECT_GT(h.weight, 0.0);
			EXPECT_GE(h.angle, 0.0);
			EXPECT_LE(h.angle, 90.0);
			EXPECT_GT(h.energy, -p.phi_wall * 0.99);
		}
	}
}

// At normal incidence E is along B, so neither changes an ion's speed
// across B: it hits with the speed across B it started with. That's
// Rayleigh at Ti, mean v_th sqrt(pi / 2), whatever the speed along B (the
// start only rejects on that), so speed * sin(angle) has that mean. (An
// ion that starts almost at rest along B can take longer than the
// sampler follows it, so a few may be lost.)
TEST(SheathImpacts, NormalIncidence)
{
	const auto p {params_at(0.0)};
	const auto r {Sheath::sample_impacts(p, Sheath::Model::brooks, 4000,
		3)};
	const int n {static_cast<int>(r.hits.size())};
	ASSERT_GE(n, 3990);
	const double mass {mi_amu * 1.661e-27};
	const double v_th {std::sqrt(1.602e-19 * 20.0 / mass)};
	double sum {};
	for (const auto& h : r.hits)
	{
		const double v {std::sqrt(2.0 * h.energy * 1.602e-19 / mass)};
		sum += v * std::sin(h.angle * pi / 180.0);
	}
	// The mean's standard error is about 0.52 v_th / sqrt(n)
	EXPECT_NEAR(sum / n, v_th * std::sqrt(pi / 2.0),
		4.0 * 0.52 * v_th / std::sqrt(n));
}

// The magnetic presheath turns the ions from along B toward the normal,
// so on average they hit nearer the normal than the field is, for both
// profiles.
TEST(SheathImpacts, TurnedTowardNormal)
{
	for (const auto m : {Sheath::Model::brooks, Sheath::Model::borodkina})
	{
		for (const double deg : {45.0, 60.0, 85.0})
		{
			const auto r {Sheath::sample_impacts(params_at(deg), m, 500, 11)};
			double w {}, mean {};
			for (const auto& h : r.hits)
			{
				w += h.weight;
				mean += h.weight * h.angle;
			}
			mean /= w;
			EXPECT_LT(mean, deg - 5.0) << static_cast<int>(m) << " " << deg;
		}
	}
}

// Each ion has its own random stream, so the same seed gives the same
// impacts (whatever the threads), and another seed different ones.
TEST(SheathImpacts, Reproducible)
{
	const auto p {params_at(60.0)};
	const auto a {Sheath::sample_impacts(p, Sheath::Model::brooks, 100, 5)};
	const auto b {Sheath::sample_impacts(p, Sheath::Model::brooks, 100, 5)};
	const auto c {Sheath::sample_impacts(p, Sheath::Model::brooks, 100, 6)};
	ASSERT_EQ(a.hits.size(), b.hits.size());
	for (std::size_t k {}; k < a.hits.size(); ++k)
	{
		EXPECT_EQ(a.hits[k].energy, b.hits[k].energy);
		EXPECT_EQ(a.hits[k].angle, b.hits[k].angle);
		EXPECT_EQ(a.hits[k].weight, b.hits[k].weight);
	}
	ASSERT_FALSE(c.hits.empty());
	EXPECT_NE(a.hits[0].energy, c.hits[0].energy);
}

// Nothing to push without a plasma, and Params carries what the sampler
// needs (Ti, B, the ion mass)
TEST(SheathImpacts, Params)
{
	const auto p {params_at(30.0, 10.0)};
	EXPECT_EQ(p.ti, 10.0);
	EXPECT_EQ(p.bmag, 2.0);
	EXPECT_EQ(p.ion_mass_amu, mi_amu);
	const auto r {Sheath::sample_impacts(Sheath::Params {},
		Sheath::Model::brooks, 10, 1)};
	EXPECT_TRUE(r.hits.empty());
	EXPECT_EQ(r.lost, 10);
}
