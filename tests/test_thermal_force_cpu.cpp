/*
* Tests of the heat-flux skew used for the thermal force: the size of the
* skew eps and the Grad distribution it defines.
*/
#include <cmath>
#include <gtest/gtest.h>

#include "collisions.h"
#include "constants.h"

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
