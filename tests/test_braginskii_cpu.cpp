/*
* Tests for the Braginskii collision frequencies and parallel heat
* conduction (braginskii.cpp): the frequencies against the NRL
* Formulary's Braginskii collision times, the Coulomb logarithms' ranges
* and floors, the conductivity's temperature scaling and the flux
* limiter.
*/

#include <cmath>
#include <gtest/gtest.h>

#include "braginskii.h"

namespace
{
	constexpr double n {1e19};		// m^-3
	constexpr double aa {2.0};		// deuterium
	constexpr double mi {aa * Braginskii::mp};
}

// nu_ei and nu_ii are 1 / tau_e and 1 / tau_i of the NRL Formulary,
// tau_e = 3.44e5 T^1.5 / (n lnL) and tau_i = 2.09e7 T^1.5 sqrt(A) / (n lnL)
// (n in cm^-3), to the Formulary's 3 digits. nu_ei has small corrections
// from the ions' speed and mass (me / mi), so its tolerance is a little
// looser than the rounding.
TEST(Braginskii, MatchesNrlCollisionTimes)
{
	for (const double t : {5.0, 50.0, 500.0})
	{
		const double lnl_ei {Braginskii::coulomb_log_ei(n, n, t, t, 1.0,
			aa)};
		const double tau_e {3.44e5 * std::pow(t, 1.5) / (n * 1e-6 * lnl_ei)};
		EXPECT_NEAR(Braginskii::nu_ei(n, n, t, t, 1.0, aa) * tau_e, 1.0,
			4e-3) << "T = " << t;

		const double lnl_ii {Braginskii::coulomb_log_ii(n, t, 1.0, aa, n, t,
			1.0, aa)};
		const double tau_i {2.09e7 * std::pow(t, 1.5) * std::sqrt(aa)
			/ (n * 1e-6 * lnl_ii)};
		EXPECT_NEAR(Braginskii::nu_ii(n, t, 1.0, aa, n, t, 1.0, aa) * tau_i,
			1.0, 4e-3) << "T = " << t;
	}
}

// Like-particle collisions are slower than electron-ion ones by sqrt(2)
// at the same Coulomb logarithm (nu_ee = nu_ei / sqrt(2) for Z = 1)
TEST(Braginskii, ElectronElectronRatio)
{
	const double t {50.0};
	const double ee {Braginskii::nu_ee(n, t)
		/ Braginskii::coulomb_log_ee(n, t)};
	const double ei {Braginskii::nu_ei(n, n, t, 0.0, 1.0, 1e6)
		/ Braginskii::coulomb_log_ei(n, n, t, 0.0, 1.0, 1e6)};
	EXPECT_NEAR(ee / ei, 1.0 / std::sqrt(2.0), 1e-6);
}

// The electron-ion Coulomb logarithm's ranges: a fixed 10 for cold or
// thin plasma, the ion form for Te < Ti me / mi, and two electron forms
// that meet at Te = e^2 Z^2
TEST(Braginskii, ElectronIonCoulombLogRanges)
{
	EXPECT_EQ(Braginskii::coulomb_log_ei(n, n, 0.05, 1.0, 1.0, aa), 10.0);
	EXPECT_EQ(Braginskii::coulomb_log_ei(1e9, n, 50.0, 1.0, 1.0, aa), 10.0);

	// Te < Ti me / mi: needs very hot ions
	const double ti {2000.0}, te {0.2};
	ASSERT_LT(te, ti * Braginskii::me / mi);
	EXPECT_NEAR(Braginskii::coulomb_log_ei(n, n, te, ti, 1.0, aa),
		23.0 - 0.5 * std::log(n) + 1.5 * std::log(ti) - std::log(aa), 1e-12);

	// Continuous where the two electron forms meet, for Z = 1 and 2
	for (const double z : {1.0, 2.0})
	{
		const double edge {std::exp(2.0) * z * z};
		EXPECT_NEAR(Braginskii::coulomb_log_ei(n, n, edge * (1.0 - 1e-9),
			1.0, z, aa), Braginskii::coulomb_log_ei(n, n, edge, 1.0, z, aa),
			1e-6);
	}
	EXPECT_NEAR(Braginskii::coulomb_log_ei(n, n, 50.0, 1.0, 1.0, aa),
		31.0 - 0.5 * std::log(n) + std::log(50.0), 1e-12);
}

// Coulomb logarithms and frequencies stay finite with no plasma, and
// the frequencies are zero with no density
TEST(Braginskii, Floors)
{
	for (const double v : {Braginskii::coulomb_log_ee(0.0, 0.0),
		Braginskii::coulomb_log_ii(0.0, 0.0, 1.0, aa, 0.0, 0.0, 1.0, aa)})
	{
		EXPECT_TRUE(std::isfinite(v));
	}
	EXPECT_EQ(Braginskii::nu_ee(0.0, 0.0), 0.0);
	EXPECT_EQ(Braginskii::nu_ei(0.0, 0.0, 0.0, 0.0, 1.0, aa), 0.0);
	EXPECT_EQ(Braginskii::kappa_par(3.9, 0.0, 50.0, 1e5, mi), 0.0);
	EXPECT_EQ(Braginskii::kappa_par(3.9, n, 50.0, 0.0, mi), 0.0);
}

// Spitzer's electron conductivity: Braginskii's 3.16 with electron-ion
// collisions alone is kappa_e = 3.16 e^2 n T tau_e / me, about
// 3.06e4 Te^2.5 / lnL W / m / eV with the NRL Formulary's tau_e (and
// 3.09e4 in Stangeby's book)
TEST(Braginskii, SpitzerConductivity)
{
	for (const double t : {10.0, 50.0, 200.0})
	{
		const double nu {Braginskii::nu_ei(n, n, t, t, 1.0, aa)};
		const double lnl {Braginskii::coulomb_log_ei(n, n, t, t, 1.0, aa)};
		EXPECT_NEAR(Braginskii::kappa_par(3.16, n, t, nu, Braginskii::me)
			/ (3.06e4 * std::pow(t, 2.5) / lnl), 1.0, 5e-3) << "T = " << t;
	}
}

// kappa = coef n e^2 T / (m nu), so with nu ~ n lnL / T^1.5 it goes as
// T^2.5 / lnL and doesn't depend on the density except through lnL
TEST(Braginskii, ConductivityScaling)
{
	auto kappa_e = [](double ne, double t)
	{
		const double nu {Braginskii::nu_ee(ne, t)
			+ Braginskii::nu_ei(ne, ne, t, t, 1.0, aa)};
		return Braginskii::kappa_par(Braginskii::kappa_coef_e_default, ne,
			t, nu, Braginskii::me);
	};
	EXPECT_NEAR(kappa_e(n, 50.0), Braginskii::kappa_coef_e_default * n
		* 50.0 * Braginskii::qe * Braginskii::qe / (Braginskii::me * (Braginskii::nu_ee(n,
		50.0) + Braginskii::nu_ei(n, n, 50.0, 50.0, 1.0, aa))),
		1e-12 * kappa_e(n, 50.0));

	const double lnl_ratio {Braginskii::coulomb_log_ei(n, n, 100.0, 100.0,
		1.0, aa) / Braginskii::coulomb_log_ei(n, n, 50.0, 50.0, 1.0, aa)};
	const double ratio {kappa_e(n, 100.0) / kappa_e(n, 50.0)};
	EXPECT_NEAR(ratio * lnl_ratio, std::pow(2.0, 2.5), 0.02);
	EXPECT_NEAR(kappa_e(1e20, 50.0) / kappa_e(n, 50.0), 1.0, 0.1);
}

// The flux limiter does nothing without alpha, keeps the heat flux under
// alpha times the free-streaming flux, and barely changes a small one
TEST(Braginskii, FluxLimiter)
{
	const double t {50.0};
	const double kappa {1e3};
	const double q_fl {0.2 * n * Braginskii::qe * t
		* std::sqrt(Braginskii::qe * t / Braginskii::me)};
	EXPECT_EQ(Braginskii::limit_kappa(kappa, 1e6, n, t, Braginskii::me,
		-1.0), kappa);
	for (const double g : {1e3, 1e6, 1e9})
	{
		const double k {Braginskii::limit_kappa(kappa, g, n, t,
			Braginskii::me, 0.2)};
		EXPECT_LT(k * g, q_fl);
		EXPECT_NEAR(k, kappa / (1.0 + kappa * g / q_fl), 1e-9 * kappa);
	}
	// q_SH / q_fl is about 2e-5 here
	EXPECT_NEAR(Braginskii::limit_kappa(kappa, 1.0, n, t, Braginskii::me,
		0.2), kappa, 1e-4 * kappa);
}
