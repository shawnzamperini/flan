/*
* Tests of the heat-flux skew used for the thermal force: the size of the
* skew eps, its shape, and sampling the distribution they define.
*/
#include <array>
#include <cmath>
#include <cstdlib>
#include <gtest/gtest.h>
#include <omp.h>
#include <string>
#include <vector>

#include "background.h"
#include "collisions.h"
#include "constants.h"
#include "flan_types.h"
#include "impurity_stats.h"
#include "read_bkg.h"
#include "slots.h"
#include "vectors.h"
#include "zones.h"
#include "options.h"
#include "pcg32.h"

namespace
{
	const Collisions::SkewShape grad_shape {Collisions::skew_shape(1.0)};
	const Collisions::SkewShape ion_shape {Collisions::skew_shape(
		Collisions::skew_push_ion)};
	const Collisions::SkewShape elec_shape {Collisions::skew_shape(
		Collisions::skew_push_elec)};

	// Moments of the skewed distribution f_M [1 + eps w_par (a + b w^2 +
	// c w^4)] in units of v_th, integrated over (w_par, w_perp) with
	// d^3w = 2 pi w_perp dw_perp dw_par and the midpoint rule. f_M is the
	// unit Maxwellian exp(-w^2 / 2) / (2 pi)^(3 / 2).
	struct Moments
	{
		double n {};		// int f
		double u {};		// int w_par f
		double T {};		// int w^2 f / 3
		double q {};		// int 1/2 w^2 w_par f
	};

	Moments skew_moments(const double eps,
		const Collisions::SkewShape& shape)
	{
		constexpr int nw {800};
		constexpr double wmax {12.0};
		constexpr double dpar {2.0 * wmax / nw};
		constexpr double dperp {wmax / nw};
		const double norm {std::pow(2.0 * Constants::pi, -1.5)};
		Moments m {};
		for (int a {}; a < nw; ++a)
		for (int b {}; b < nw; ++b)
		{
			const double w_par {-wmax + (a + 0.5) * dpar};
			const double w_perp {(b + 0.5) * dperp};
			const double w_sq {w_par * w_par + w_perp * w_perp};
			const double f {norm * std::exp(-0.5 * w_sq)
				* Collisions::skew_factor(eps, w_par, w_sq, shape)
				* 2.0 * Constants::pi * w_perp * dpar * dperp};
			m.n += f;
			m.u += w_par * f;
			m.T += w_sq * f / 3.0;
			m.q += 0.5 * w_sq * w_par * f;
		}
		return m;
	}

	// Average of g(w) w_par^2 (a + b w^2 + c w^4) over the unit
	// Maxwellian, where g depends on the speed only: by symmetry it is
	// 1/3 of the average of g(w) w^2 (a + ...), a radial integral with
	// the speed distribution sqrt(2 / pi) w^2 exp(-w^2 / 2)
	template <typename G>
	double shape_average(const Collisions::SkewShape& shape, G g)
	{
		constexpr int nw {200000};
		constexpr double wmax {14.0};
		constexpr double dw {wmax / nw};
		double sum {};
		for (int i {}; i < nw; ++i)
		{
			const double w {(i + 0.5) * dw};
			const double w2 {w * w};
			sum += g(w) * w2 * (shape.a + shape.b * w2 + shape.c * w2 * w2)
				* std::sqrt(2.0 / Constants::pi) * w2 * std::exp(-0.5 * w2)
				* dw;
		}
		return sum / 3.0;
	}
}

// eps = q_par / (n T v_th), with v_th = sqrt(T / m)
TEST(ThermalForce, SkewEps)
{
	const double n {1e19};
	const double T {10.0};
	const double m {2.014 * Constants::amu_to_kg};
	const double T_J {T * Constants::ev_to_j};
	const double v_th {std::sqrt(T_J / m)};
	const double q {0.07 * n * T_J * v_th};
	EXPECT_NEAR(Collisions::skew_eps(q, n, T, m), 0.07, 1e-12);
	EXPECT_NEAR(Collisions::skew_eps(-q, n, T, m), -0.07, 1e-12);
	EXPECT_EQ(Collisions::skew_eps(0.0, n, T, m), 0.0);
}

// The shape's three conditions, per unit eps: no flow, a heat flux of 1,
// and a friction moment <w_par^2 (...) / w^3> (the push on a slow, heavy
// impurity) of push times the Grad shape's, -0.2 sqrt(2 / pi). push = 1
// is the Grad shape w_par (w^2 / 5 - 1) itself.
TEST(ThermalForce, SkewShape)
{
	EXPECT_DOUBLE_EQ(grad_shape.a, -1.0);
	EXPECT_DOUBLE_EQ(grad_shape.b, 0.2);
	EXPECT_DOUBLE_EQ(grad_shape.c, 0.0);
	EXPECT_NEAR(grad_shape.bound, 20.0, 1e-12);

	// The closure values: 2.65 against Grad's 3.31 for ions, and 0.71
	// against 1.90 for electrons
	EXPECT_NEAR(Collisions::skew_push_ion, 0.800, 1e-3);
	EXPECT_NEAR(Collisions::skew_push_elec, 0.3745, 1e-4);

	const double grad_friction {-0.2 * std::sqrt(2.0 / Constants::pi)};
	for (const double push : {1.0, Collisions::skew_push_ion,
		Collisions::skew_push_elec, 0.0, 1.5})
	{
		const auto shape {Collisions::skew_shape(push)};
		EXPECT_NEAR(shape_average(shape, [](double) {return 1.0;}), 0.0,
			1e-8) << push;
		EXPECT_NEAR(shape_average(shape, [](double w)
			{return 0.5 * w * w;}), 1.0, 1e-8) << push;
		EXPECT_NEAR(shape_average(shape, [](double w)
			{return 1.0 / (w * w * w);}), push * grad_friction, 1e-8)
			<< push;
	}
}

// The skew keeps the density, flow and temperature of the Maxwellian and
// adds a heat flux of exactly eps n T v_th (so eps from skew_eps gives
// back q_par), whatever the shape. Integrating the full factor, negative
// parts included, so this holds for any eps. The tolerance is the
// midpoint rule's error.
TEST(ThermalForce, SkewMoments)
{
	for (const auto* shape : {&grad_shape, &ion_shape, &elec_shape})
	for (const double eps : {0.0, 0.07, -0.3, 1.0})
	{
		const auto m {skew_moments(eps, *shape)};
		EXPECT_NEAR(m.n, 1.0, 1e-4) << eps;
		EXPECT_NEAR(m.u, 0.0, 1e-4) << eps;
		EXPECT_NEAR(m.T, 1.0, 1e-4) << eps;
		EXPECT_NEAR(m.q, eps, 2e-4) << eps;
	}
}

// Fast particles lean along the heat flux and slow ones against it. For
// the Grad shape the factor changes sign at w^2 = 5.
TEST(ThermalForce, SkewFactor)
{
	const double eps {0.1};
	for (const auto* shape : {&grad_shape, &ion_shape, &elec_shape})
	{
		EXPECT_GT(Collisions::skew_factor(eps, 4.0, 16.0, *shape), 1.0);
		EXPECT_LT(Collisions::skew_factor(eps, -4.0, 16.0, *shape), 1.0);
		EXPECT_LT(Collisions::skew_factor(eps, 1.0, 1.0, *shape), 1.0);
		EXPECT_GT(Collisions::skew_factor(eps, -1.0, 1.0, *shape), 1.0);
		EXPECT_DOUBLE_EQ(Collisions::skew_factor(0.0, 3.0, 9.0, *shape),
			1.0);
	}
	EXPECT_DOUBLE_EQ(Collisions::skew_factor(eps, 2.0, 5.0, grad_shape),
		1.0);

	// The slowest electrons lean much less than the Grad shape's: the
	// electron shape's lean sits at higher speeds
	EXPECT_LT(std::abs(elec_shape.a), 0.2 * std::abs(grad_shape.a));
}

namespace
{
	// Moments of velocities from sample_bkg_velocity, in units of
	// v_th relative to the flow: mean w, temperature (mean w^2 / 3), and the
	// heat flux along b and across it, mean 1/2 w^2 w (which is
	// q / (n T v_th), i.e., eps if the skew is reproduced)
	struct Sampled
	{
		double w_par {};
		double w_perp {};
		double T {};
		double q_par {};
		double q_perp {};
		Collisions::SkewCounts counts {};
	};

	Sampled sample_moments(const double eps,
		const Collisions::SkewShape& shape, const int n)
	{
		const double T {20.0};
		const double m {2.014 * Constants::amu_to_kg};
		const double v_th {std::sqrt(T * Constants::ev_to_j / m)};
		const double u[3] {3e4, -1e4, 2e4};
		const double b[3] {1.0 / 3.0, 2.0 / 3.0, -2.0 / 3.0};
		const double e_perp[3] {2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0};	// . b = 0
		pcg32 rng(42, 7);
		Sampled s {};
		for (int k {}; k < n; ++k)
		{
			const auto [vX, vY, vZ] = Collisions::sample_bkg_velocity(
				T, u[0], u[1], u[2], m, eps, shape, b[0], b[1], b[2], rng,
				s.counts);
			const double w[3] {(vX - u[0]) / v_th, (vY - u[1]) / v_th,
				(vZ - u[2]) / v_th};
			const double w_par {w[0] * b[0] + w[1] * b[1] + w[2] * b[2]};
			const double w_perp {w[0] * e_perp[0] + w[1] * e_perp[1]
				+ w[2] * e_perp[2]};
			const double w_sq {w[0] * w[0] + w[1] * w[1] + w[2] * w[2]};
			s.w_par += w_par / n;
			s.w_perp += w_perp / n;
			s.T += w_sq / 3.0 / n;
			s.q_par += 0.5 * w_sq * w_par / n;
			s.q_perp += 0.5 * w_sq * w_perp / n;
		}
		return s;
	}
}

// The bound on the skew factor holds for |w| <= 5, and is reached at
// w_par = w = 5 for these shapes
TEST(ThermalForce, SkewBound)
{
	EXPECT_DOUBLE_EQ(Collisions::skew_bound(0.07, grad_shape),
		1.0 + 20.0 * 0.07);
	for (const auto* shape : {&grad_shape, &ion_shape, &elec_shape})
	for (const double eps : {0.07, -0.3, 1.0})
	{
		const double M {Collisions::skew_bound(eps, *shape)};
		const double w {eps > 0.0 ? 5.0 : -5.0};
		EXPECT_NEAR(Collisions::skew_factor(eps, w, 25.0, *shape), M,
			1e-12 * M);
		for (double wp {-5.0}; wp <= 5.0; wp += 0.25)
		for (double wq {0.0}; wp * wp + wq * wq <= 25.0; wq += 0.25)
		{
			EXPECT_LE(Collisions::skew_factor(eps, wp, wp * wp + wq * wq,
				*shape), M * (1.0 + 1e-12));
		}
	}
}

// Sampled background particles keep the flow and temperature, and carry
// the requested heat flux along B and none across it. Clipping the
// negative tail changes the moments a little when the skew is strong.
// With the ion shape: nothing measurable at eps = 0.07, but at |eps| = 0.3
// (0.28% of draws clipped) the heat flux is 0.86 of eps, the mean shifts
// by 0.004 v_th toward the hot side (against q) and T rises by 0.7%. With
// the electron shape at eps = 0.1 (0.045% clipped), the heat flux is 0.94
// of eps. The statistical errors are about 5e-4 on the mean and 0.0015 on
// the heat flux.
TEST(ThermalForce, SkewedSampler)
{
	const int n {4000000};
	struct Case
	{
		const Collisions::SkewShape* shape;
		double eps, q_lo, q_hi, w_par, w_tol, T, T_tol, clip_lo, clip_hi;
	};
	for (const auto c : std::vector<Case> {
		{&ion_shape, 0.07, 0.066, 0.073, 0.0, 2e-3, 1.0, 3e-3, 0.0, 1e-4},
		{&ion_shape, -0.3, -0.264, -0.250, 0.0040, 2e-3, 1.0071, 3e-3,
			2.3e-3, 3.3e-3},
		{&elec_shape, 0.02, 0.016, 0.024, 0.0, 2e-3, 1.0, 3e-3, 0.0, 1e-4},
		{&elec_shape, 0.1, 0.090, 0.099, 0.0, 2e-3, 1.0, 3e-3, 3e-4,
			6e-4}})
	{
		const double eps {c.eps};
		const auto s {sample_moments(eps, *c.shape, n)};
		EXPECT_NEAR(s.w_par, c.w_par, c.w_tol) << eps;
		EXPECT_NEAR(s.w_perp, 0.0, 2e-3) << eps;
		EXPECT_NEAR(s.T, c.T, c.T_tol) << eps;
		EXPECT_GT(s.q_par, c.q_lo) << eps;
		EXPECT_LT(s.q_par, c.q_hi) << eps;
		EXPECT_NEAR(s.q_perp, 0.0, 6e-3) << eps;

		// About M draws per sample, and the share of draws clipped
		EXPECT_EQ(s.counts.samples, n);
		const double M {Collisions::skew_bound(eps, *c.shape)};
		EXPECT_NEAR(double(s.counts.draws) / n, M, 0.02 * M) << eps;
		const double clipped {double(s.counts.negative) / s.counts.draws};
		EXPECT_GE(clipped, c.clip_lo) << eps;
		EXPECT_LE(clipped, c.clip_hi) << eps;
		EXPECT_EQ(s.counts.eps_capped, 0);
	}

	// |eps| > 1 is capped at 1 and counted
	const auto s {sample_moments(3.0, ion_shape, 1000)};
	EXPECT_EQ(s.counts.eps_capped, 1000);
	const double M {Collisions::skew_bound(1.0, ion_shape)};
	EXPECT_NEAR(double(s.counts.draws) / 1000, M, 0.15 * M);
}

// Without a skew, the sampler is the plain drifting Maxwellian, with the
// same random numbers as before the skew was added (each component the
// flow plus rng.normal times v_th), and nothing is counted.
//
// The velocities are compared to within rounding, not bit for bit: the
// library may be compiled with -march=native (Release builds) and this
// test without it, and then u + v_th * w can be a fused multiply-add on
// one side and not the other, which changes the last bits. Whether the
// same random numbers were used is checked exactly instead, by drawing
// from both generators afterward.
TEST(ThermalForce, UnskewedSampler)
{
	const double T {20.0};
	const double m {2.014 * Constants::amu_to_kg};
	const double v_th {std::sqrt(T * Constants::ev_to_j / m)};
	const double u[3] {3e4, -1e4, 2e4};
	pcg32 rng(3, 5);
	pcg32 old_rng(3, 5);
	Collisions::SkewCounts counts {};
	for (int k {}; k < 1000; ++k)
	{
		const auto [vX, vY, vZ] = Collisions::sample_bkg_velocity(T, u[0],
			u[1], u[2], m, 0.0, ion_shape, 0.0, 0.0, 1.0, rng, counts);
		const double v[3] {vX, vY, vZ};
		for (int c {}; c < 3; ++c)
		{
			EXPECT_NEAR(v[c], old_rng.normal(u[c], v_th),
				1e-13 * (std::abs(u[c]) + v_th)) << k << " " << c;
		}
	}
	EXPECT_EQ(rng.next_uint(), old_rng.next_uint());
	EXPECT_EQ(counts.samples, 0);
	EXPECT_EQ(counts.draws, 0);
}

// The thermal force is on by default. Without heat fluxes it only warns
// (any device); with only one heat flux, or with heat fluxes on the GPU, it
// stops the run; off checks nothing.
TEST(ThermalForce, Checks)
{
	Background::Background bkg {};
	Options::Options opts {};
	EXPECT_EQ(opts.imp_thermal_force_int(), 1);

	testing::internal::CaptureStdout();
	Collisions::check_thermal_force(bkg, opts);
	EXPECT_NE(testing::internal::GetCapturedStdout().find("no parallel heat "
		"fluxes, so the thermal force isn't"), std::string::npos);
	opts.set_use_gpu("cuda");
	testing::internal::CaptureStdout();
	Collisions::check_thermal_force(bkg, opts);
	EXPECT_NE(testing::internal::GetCapturedStdout().find("Warning!"),
		std::string::npos);

	// One zone, one cell and frame, so the background can hold heat fluxes
	Zones::Zones zones(0);
	zones.add_zone({0.5}, {0.5}, {0.5}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0});
	zones.set_face(0, Zones::Face::xmin, Zones::face_absorbing);
	bkg.move_into_zones(zones);
	auto heat_flux = [&]()
	{
		return Vectors::ZonedVector<BkgFPType>(1, bkg.get_zones().get_nx_all(),
			bkg.get_zones().get_ny_all(), bkg.get_zones().get_nz_all(), 0);
	};
	auto q_i {heat_flux()};
	bkg.move_into_q_par_i(q_i);
	EXPECT_EXIT(Collisions::check_thermal_force(bkg, opts),
		testing::ExitedWithCode(EXIT_FAILURE), "only has q_par_i");
	auto q_e {heat_flux()};
	bkg.move_into_q_par_e(q_e);
	EXPECT_EXIT(Collisions::check_thermal_force(bkg, opts),
		testing::ExitedWithCode(EXIT_FAILURE), "only runs on the CPU");

	opts.set_use_gpu("off");
	testing::internal::CaptureStdout();
	Collisions::check_thermal_force(bkg, opts);
	EXPECT_EQ(testing::internal::GetCapturedStdout(), "");

	opts.set_imp_thermal_force("off");
	opts.set_use_gpu("cuda");
	Collisions::check_thermal_force(bkg, opts);
}

namespace
{
	// Tungsten at rest in the thermal_force test background (uniform
	// ne = 1e20 m^-3, Te = Ti = 1 eV, B along +X, heat fluxes in +X with
	// eps = 0.1 for both species), and what the collision model does to it
	struct ForceCase
	{
		Options::Options opts {};
		Background::Background bkg {};
		static constexpr double imp_mass {183.84 * Constants::amu_to_kg};
		static constexpr double ne {1e20};
		static constexpr double T {1.0};

		ForceCase()
		{
			opts.set_bkg_source("test");
			opts.set_test_opt("thermal_force");
			bkg = Background::read_bkg(opts);
		}

		// Mass (kg) of the background species collided with
		double bkg_mass(const bool elec) const
		{
			return (elec ? bkg.get_elec_mass_amu() : bkg.get_ion_mass_amu())
				* Constants::amu_to_kg;
		}

		// Thermal speed sqrt(T / m) of that species (m/s)
		double v_th(const bool elec) const
		{
			return std::sqrt(T * Constants::ev_to_j / bkg_mass(elec));
		}

		// Skew of that species at the middle of the slab
		double eps(const bool elec) const
		{
			const double t {bkg.get_times()[0]};
			const double q_par {elec
				? bkg.interp_q_par_e(0, t, 0.025, 0.0, 0.0)
				: bkg.interp_q_par_i(0, t, 0.025, 0.0, 0.0)};
			return Collisions::skew_eps(q_par, ne, T, bkg_mass(elec));
		}

		// Nanbu's s (how collisional a time step is) per second, for a
		// background particle at the thermal speed hitting an impurity of
		// charge Z at rest: s = s_rate dt (v_th / g)^3 at relative speed
		// g. The same formulas as nanbu_calc_s, written out again here.
		double s_rate(const bool elec, const int Z) const
		{
			const double e {Constants::charge_e};
			const double mb {bkg_mass(elec)};
			const double mu {imp_mass * mb / (imp_mass + mb)};
			const double T_J {T * Constants::ev_to_j};
			const double b0 {Z * e * e
				/ (2.0 * Constants::pi * Constants::eps0 * 3.0 * T_J)};
			const double debye {std::sqrt(Constants::eps0 * T_J
				/ (ne * e * e))};
			const double ln_alpha {std::max(std::log(debye / b0), 1.0)};
			const double sq {Z * e * e / (Constants::eps0 * mu)};
			const double vt {v_th(elec)};
			return ln_alpha / (4.0 * Constants::pi) * sq * sq * ne
				/ (vt * vt * vt);
		}
	};

	// Mean and standard error of the change in velocity along B (+X) and
	// across it (+Y) that one call of nanbu_coll gives impurities at
	// rest, in m/s, from n_slots * n_calls collisions
	struct Kick
	{
		double par {};
		double par_err {};
		double perp {};
		double perp_err {};
	};

	Kick measure_kick(const ForceCase& fc, const bool elec, const int Z,
		const double dt, const int n_slots, const int n_calls)
	{
		Slots::Slots slots(n_slots);
		slots.set_mass(ForceCase::imp_mass);
		for (int i {}; i < n_slots; ++i)
		{
			slots.set_zone(i, 0);
			slots.set_t(i, fc.bkg.get_times()[0]);
			slots.set_x(i, 0.025);
			slots.set_y(i, 0.0);
			slots.set_z(i, 0.0);
			slots.set_q(i, Z);
			slots.set_state(i, 0);
		}

		Impurity::Statistics imp_stats {};
		std::vector<pcg32> rngs {};
		for (int i {}; i < omp_get_max_threads(); ++i)
			rngs.emplace_back(1234, 2 * i + (elec ? 1 : 0));

		double sum_par {}, sum_par_sq {}, sum_perp {}, sum_perp_sq {};
		for (int c {}; c < n_calls; ++c)
		{
			for (int i {}; i < n_slots; ++i)
			{
				slots.set_vX(i, 0.0);
				slots.set_vY(i, 0.0);
				slots.set_vZ(i, 0.0);
			}
			Collisions::nanbu_coll(slots, fc.bkg, fc.opts, elec, dt,
				imp_stats, rngs);
			for (int i {}; i < n_slots; ++i)
			{
				const double vX {slots.vX()[i]};
				const double vY {slots.vY()[i]};
				sum_par += vX;
				sum_par_sq += vX * vX;
				sum_perp += vY;
				sum_perp_sq += vY * vY;
			}
		}

		const double N {static_cast<double>(n_slots) * n_calls};
		Kick k {};
		k.par = sum_par / N;
		k.par_err = std::sqrt((sum_par_sq / N - k.par * k.par) / N);
		k.perp = sum_perp / N;
		k.perp_err = std::sqrt((sum_perp_sq / N - k.perp * k.perp) / N);
		return k;
	}
}

// The collision model gives the thermal force: with a heat flux, an
// impurity at rest gets a net kick along B, opposite to the heat flux
// (toward where the plasma is hotter), of the size the model should give.
//
// In one Nanbu step a background particle with velocity v_b turns the
// relative velocity by an angle chi, which changes the velocity of an
// impurity at rest by frac v_b (1 - cos chi) on average over the random
// azimuth, with frac = m_b / (m_imp + m_b). Nanbu's model has <cos chi> =
// exp(-s), and s = s0 / w^3 for a background particle at w thermal speeds
// (s0 = s_rate dt). Averaging over the skewed distribution, the Maxwellian
// part cancels by symmetry and
//
//   <dv_par> = frac v_th eps <w_par^2 (a + b w^2 + c w^4)
//                             (1 - exp(-s0 / w^3))>,
//
// exact at any time step and linear in eps. Slow particles (large s) turn
// all the way and fast ones barely, and since the slow ones lean against
// the heat flux, the kick is against it. Checked for both species; the
// last check is that it is a property of the heat flux (no heat flux, no
// kick).
//
// The kick is small next to the random part of a collision, so resolving it
// takes millions of collisions. An optimized build does 2e7 for each
// species at two time steps (s0 = 0.3 and 3), which resolves the kick to
// 2 - 5%. A build without optimization is about 50 times slower, so it
// does 2e6 at s0 = 3 only (5 - 7%), enough to catch a wrong sign or a
// badly wrong size.
TEST(ThermalForce, NanbuKick)
{
	const ForceCase fc {};
	constexpr int Z {2};
	constexpr int n_slots {200000};
#ifdef __OPTIMIZE__
	constexpr int n_calls {100};
	const std::vector<double> s0s {0.3, 3.0};
#else
	constexpr int n_calls {10};
	const std::vector<double> s0s {3.0};
#endif

	for (const bool elec : {false, true})
	{
		const Collisions::SkewShape& shape {elec ? elec_shape : ion_shape};
		const double eps {fc.eps(elec)};
		EXPECT_NEAR(eps, 0.1, 1e-12);
		const double mb {fc.bkg_mass(elec)};
		const double frac {mb / (ForceCase::imp_mass + mb)};

		for (const double s0 : s0s)
		{
			const double dt {s0 / fc.s_rate(elec, Z)};
			const double expected {frac * fc.v_th(elec) * eps
				* shape_average(shape, [s0](const double w)
				{return 1.0 - std::exp(-s0 / (w * w * w));})};
			const Kick k {measure_kick(fc, elec, Z, dt, n_slots, n_calls)};

			SCOPED_TRACE((elec ? "electrons, s0 = " : "ions, s0 = ")
				+ std::to_string(s0));
			std::cout << (elec ? "electrons" : "ions") << ", s0 = " << s0
				<< ": dv_par = " << k.par << " +/- " << k.par_err
				<< " m/s, expected " << expected << '\n';

			// Against the heat flux (+X), and resolved well enough for
			// the comparison to mean something
			EXPECT_LT(expected, 0.0);
			EXPECT_LT(k.par, 0.0);
			EXPECT_LT(k.par_err, 0.1 * std::abs(expected));

			// Agrees to within 4 standard errors, plus 3% for what the
			// sampler clips at eps = 0.1 and Nanbu's table of A(s)
			EXPECT_NEAR(k.par, expected,
				4.0 * k.par_err + 0.03 * std::abs(expected));

			// Nothing across B
			EXPECT_NEAR(k.perp, 0.0, 4.0 * k.perp_err);
		}
	}

	// Option off: Maxwellian background, no kick
	ForceCase off {};
	off.opts.set_imp_thermal_force("off");
	const double dt {1.0 / off.s_rate(false, Z)};
	const Kick k {measure_kick(off, false, Z, dt, n_slots, n_calls / 5)};
	EXPECT_NEAR(k.par, 0.0, 4.0 * k.par_err);
}

// In the limit of a small time step the kick is a force, F = m_imp
// <dv_par> / dt, and with the heat flux of a collisional plasma, q_par =
// -kappa grad T, it is the textbook thermal force: 2.65 Z^2 grad T from
// the ions (for a heavy impurity; a factor m_i / mu from the reduced mass
// mu otherwise) and 0.71 Z^2 grad T from the electrons. This ties
// together the Coulomb rate in Nanbu's s, Braginskii's conductivities
// (kappa = 3.906 n T tau_i / m_i and 3.16 n T tau_e / m_e) and the fitted
// shapes of the skew.
TEST(ThermalForce, SmallStepLimit)
{
	const ForceCase fc {};
	constexpr int Z {3};
	constexpr double grad_T {-50.0 * Constants::ev_to_j};	// J/m, hot at -X
	const double e {Constants::charge_e};
	const double T_J {ForceCase::T * Constants::ev_to_j};

	for (const bool elec : {false, true})
	{
		const Collisions::SkewShape& shape {elec ? elec_shape : ion_shape};
		const double mb {fc.bkg_mass(elec)};
		const double mu {ForceCase::imp_mass * mb / (ForceCase::imp_mass + mb)};

		// Coulomb logarithm as in s_rate, so that it cancels
		const double b0 {Z * e * e
			/ (2.0 * Constants::pi * Constants::eps0 * 3.0 * T_J)};
		const double debye {std::sqrt(Constants::eps0 * T_J
			/ (ForceCase::ne * e * e))};
		const double ln_alpha {std::max(std::log(debye / b0), 1.0)};

		// Braginskii's collision times (s) and conductivities
		// (1/(m s)): tau_i has a sqrt(2) that tau_e doesn't
		const double tau {(elec ? 6.0 * std::sqrt(2.0) : 12.0)
			* std::pow(Constants::pi, 1.5) * Constants::eps0
			* Constants::eps0 * std::sqrt(mb) * std::pow(T_J, 1.5)
			/ (ln_alpha * e * e * e * e * ForceCase::ne)};
		const double kappa {(elec ? 3.16 : 3.906) * ForceCase::ne * T_J
			* tau / mb};
		const double q_par {-kappa * grad_T};
		const double eps {Collisions::skew_eps(q_par, ForceCase::ne,
			ForceCase::T, mb)};

		// Kick per second in the small-step limit, 1 - exp(-s) -> s
		const double dv_dt {mb / (ForceCase::imp_mass + mb) * fc.v_th(elec)
			* eps * fc.s_rate(elec, Z) * shape_average(shape,
			[](const double w) {return 1.0 / (w * w * w);})};
		const double force {ForceCase::imp_mass * dv_dt};

		const double coeff {elec ? 0.71 : 2.6517 * mb / mu};
		EXPECT_NEAR(force / (Z * Z * grad_T), coeff, 2e-3 * coeff)
			<< (elec ? "electrons" : "ions");
	}
}
