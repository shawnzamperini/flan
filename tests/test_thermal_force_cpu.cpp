/*
* Tests of the heat-flux skew used for the thermal force: the size of the
* skew eps and the Grad distribution it defines.
*/
#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "collisions.h"
#include "constants.h"
#include "pcg32.h"

namespace
{
	// Moments of the Grad distribution f_M [1 + eps w_par (w^2 / 5 - 1)]
	// in units of v_th, integrated over (w_par, w_perp) with
	// d^3w = 2 pi w_perp dw_perp dw_par and the midpoint rule. f_M is the
	// unit Maxwellian exp(-w^2 / 2) / (2 pi)^(3 / 2).
	struct Moments
	{
		double n {};		// int f
		double u {};		// int w_par f
		double T {};		// int w^2 f / 3
		double q {};		// int 1/2 w^2 w_par f
	};

	Moments grad_moments(const double eps)
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
				* Collisions::grad_factor(eps, w_par, w_sq)
				* 2.0 * Constants::pi * w_perp * dpar * dperp};
			m.n += f;
			m.u += w_par * f;
			m.T += w_sq * f / 3.0;
			m.q += 0.5 * w_sq * w_par * f;
		}
		return m;
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
	EXPECT_NEAR(Collisions::grad_skew_eps(q, n, T, m), 0.07, 1e-12);
	EXPECT_NEAR(Collisions::grad_skew_eps(-q, n, T, m), -0.07, 1e-12);
	EXPECT_EQ(Collisions::grad_skew_eps(0.0, n, T, m), 0.0);
}

// The skew keeps the density, flow and temperature of the Maxwellian and
// adds a heat flux of exactly eps n T v_th (so eps from grad_skew_eps
// gives back q_par). Integrating the full factor, negative parts
// included, so this holds for any eps. The tolerance is the midpoint
// rule's error, about 1e-5.
TEST(ThermalForce, GradMoments)
{
	for (const double eps : {0.0, 0.07, -0.3, 1.0})
	{
		const auto m {grad_moments(eps)};
		EXPECT_NEAR(m.n, 1.0, 1e-4) << eps;
		EXPECT_NEAR(m.u, 0.0, 1e-4) << eps;
		EXPECT_NEAR(m.T, 1.0, 1e-4) << eps;
		EXPECT_NEAR(m.q, eps, 1e-4) << eps;
	}
}

// Fast particles lean along the heat flux and slow ones against it, and
// the factor changes sign at w^2 = 5
TEST(ThermalForce, GradFactor)
{
	const double eps {0.1};
	EXPECT_GT(Collisions::grad_factor(eps, 3.0, 9.0), 1.0);
	EXPECT_LT(Collisions::grad_factor(eps, -3.0, 9.0), 1.0);
	EXPECT_LT(Collisions::grad_factor(eps, 1.0, 1.0), 1.0);
	EXPECT_GT(Collisions::grad_factor(eps, -1.0, 1.0), 1.0);
	EXPECT_DOUBLE_EQ(Collisions::grad_factor(eps, 2.0, 5.0), 1.0);
	EXPECT_DOUBLE_EQ(Collisions::grad_factor(0.0, 3.0, 9.0), 1.0);
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

	Sampled sample_moments(const double eps, const int n)
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
				T, u[0], u[1], u[2], m, eps, b[0], b[1], b[2], rng,
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

// The bound on the Grad factor holds for |w| <= 5, and is reached there
TEST(ThermalForce, SkewBound)
{
	for (const double eps : {0.07, -0.3, 1.0})
	{
		const double M {Collisions::skew_bound(eps)};
		EXPECT_DOUBLE_EQ(M, 1.0 + 20.0 * std::abs(eps));
		const double w {eps > 0.0 ? 5.0 : -5.0};
		EXPECT_DOUBLE_EQ(Collisions::grad_factor(eps, w, 25.0), M);
		for (double wp {-5.0}; wp <= 5.0; wp += 0.25)
		for (double wq {0.0}; wp * wp + wq * wq <= 25.0; wq += 0.25)
		{
			EXPECT_LE(Collisions::grad_factor(eps, wp, wp * wp + wq * wq),
				M + 1e-12);
		}
	}
}

// Sampled background ions keep the flow and temperature, and carry the
// requested heat flux along B and none across it. Clipping the negative
// tail changes the moments a little: at eps = 0.07 nothing measurable
// (heat flux 0.99 of eps), but at |eps| = 0.3 (0.25% of the Maxwellian
// clipped) the heat flux is 0.89-0.90 of eps, the mean shifts by about
// 0.004 v_th toward the hot side (against q) and T rises by 0.5%. The
// statistical errors are about 5e-4 on the mean and 0.0015 on the heat
// flux.
TEST(ThermalForce, SkewedSampler)
{
	const int n {4000000};
	struct Case {double eps, q_lo, q_hi, w_par, w_tol, T, T_tol;};
	for (const auto c : std::vector<Case> {
		{0.07, 0.066, 0.073, 0.0, 2e-3, 1.0, 3e-3},
		{-0.3, -0.275, -0.262, 0.0039, 2e-3, 1.0053, 3e-3}})
	{
		const double eps {c.eps};
		const auto s {sample_moments(eps, n)};
		EXPECT_NEAR(s.w_par, c.w_par, c.w_tol) << eps;
		EXPECT_NEAR(s.w_perp, 0.0, 2e-3) << eps;
		EXPECT_NEAR(s.T, c.T, c.T_tol) << eps;
		EXPECT_GT(s.q_par, c.q_lo) << eps;
		EXPECT_LT(s.q_par, c.q_hi) << eps;
		EXPECT_NEAR(s.q_perp, 0.0, 6e-3) << eps;

		// About M draws per sample, and the clipped share as estimated
		// before (6e-5 at eps = 0.1, 2.5e-3 at 0.3)
		EXPECT_EQ(s.counts.samples, n);
		const double per_sample {double(s.counts.draws) / n};
		EXPECT_NEAR(per_sample, Collisions::skew_bound(eps),
			0.02 * Collisions::skew_bound(eps)) << eps;
		const double clipped {double(s.counts.negative) / s.counts.draws};
		if (std::abs(eps) < 0.1) EXPECT_LT(clipped, 1e-4);
		else EXPECT_NEAR(clipped, 2.5e-3, 0.5e-3);
		EXPECT_EQ(s.counts.eps_capped, 0);
	}

	// |eps| > 1 is capped at 1 and counted
	const auto s {sample_moments(3.0, 1000)};
	EXPECT_EQ(s.counts.eps_capped, 1000);
	EXPECT_NEAR(double(s.counts.draws) / 1000, 21.0, 2.0);
}

// Without a skew, the sampler is the plain drifting Maxwellian, with the
// same random numbers as before the skew was added (each component the
// flow plus rng.normal times v_th), and nothing is counted
TEST(ThermalForce, UnskewedSampler)
{
	const double T {20.0};
	const double m {2.014 * Constants::amu_to_kg};
	const double v_th {std::sqrt(T * Constants::ev_to_j / m)};
	pcg32 rng(3, 5);
	pcg32 old_rng(3, 5);
	Collisions::SkewCounts counts {};
	for (int k {}; k < 1000; ++k)
	{
		const auto [vX, vY, vZ] = Collisions::sample_bkg_velocity(T, 3e4,
			-1e4, 2e4, m, 0.0, 0.0, 0.0, 1.0, rng, counts);
		EXPECT_EQ(vX, old_rng.normal(3e4, v_th));
		EXPECT_EQ(vY, old_rng.normal(-1e4, v_th));
		EXPECT_EQ(vZ, old_rng.normal(2e4, v_th));
	}
	EXPECT_EQ(counts.samples, 0);
	EXPECT_EQ(counts.draws, 0);
}
