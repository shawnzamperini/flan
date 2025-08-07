/**
* @file collisions.cpp
*/
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

#include "background.h"
#include "constants.h"
#include "impurity.h"
#include "nanbu_s_a.h"
#include "options.h"
#include "random.h"
#include "utilities.h"
#include "variance_reduction.h"

namespace Collisions
{
	double calc_debye_length(const double te, const double ne)
	{
		return std::sqrt(Constants::eps0 * te * Constants::ev_to_j / 
			(ne * Constants::charge_e * Constants::charge_e)); 
	}

	std::tuple<double, double, double> sample_bkg_velocity(
		const Background::Background& bkg, const int tidx, const int xidx, 
		const int yidx, const int zidx, const double T, const double m)
	{
		// Get mean parallel flow, which is the drifting term in a drifting
		// Maxwellian. This is the ion velocity, which is assumed to be the
		// plasma velocity (so same for both ions and electrons). 
		double uX {bkg.get_uX()(tidx, xidx, yidx, zidx)};
		double uY {bkg.get_uY()(tidx, xidx, yidx, zidx)};
		double uZ {bkg.get_uZ()(tidx, xidx, yidx, zidx)};

		// Sampling from a drifting Maxwellian, which can be done by just
		// sampling from a normally distributed number with mean = drift and
		// mu = sqrt(kT/m). Unfortunately throwing away the second random
		// number here. T [eV], m [kg]
		double mu {std::sqrt(T * Constants::ev_to_j / m)};  // m/s
		double vX {};
		double vY {};
		double vZ {};
		std::tie(vX, std::ignore) = Random::get_two_norm(uX, mu);
		std::tie(vY, std::ignore) = Random::get_two_norm(uY, mu);
		std::tie(vZ, std::ignore) = Random::get_two_norm(uZ, mu);

		// Return as tuple
		return std::make_tuple(vX, vY, vZ);
	}

	// Calculate s term in Nanbu collision model
	std::tuple<double, double, double, double> nanbu_calc_s(
		Impurity::Impurity& imp, const Background::Background& bkg, 
		const int tidx, const int xidx, const int yidx, const int zidx, 
		const double imp_time_step, const double T, const double mass_kg, 
		const double Te, const double ne)
	{
		// Take random sample of background species being collided with
		auto [bkg_vX, bkg_vY, bkg_vZ] = sample_bkg_velocity(bkg, tidx, xidx, 
			yidx, zidx, T, mass_kg);

		// Relative velocity, g. Limit to a sufficiently small number to avoid
		// overflows that can occur in std::exp and std::sinh.
		double gX {imp.get_vX() - bkg_vX};
		double gY {imp.get_vY() - bkg_vY};
		double gZ {imp.get_vZ() - bkg_vZ};
		double g {std::sqrt(gX*gX + gY*gY + gZ*gZ)};
		if (g < Constants::small) g = 0.001;

		// Reduced mass of impurity and species colliding with
		double mu_ab {imp.get_mass() * mass_kg / (imp.get_mass() + mass_kg)};

		// Expectation value of g^2. Limit to a sufficiently small number to
		// avoid overflows that can occur in std::exp and std::sinh.
		double expect_g_sq {3.0 * T * Constants::ev_to_j / mu_ab};
		if (expect_g_sq < Constants::small) expect_g_sq = 0.001;

		// Expectation value of b0. Assuming singly charged background species.
		double expect_b0 {imp.get_charge() * Constants::charge_e 
			* Constants::charge_e / (2.0 * Constants::pi * Constants::eps0 
			* mu_ab * expect_g_sq)};

		// Debye length and Coulumb logarithm
		double debye_length {calc_debye_length(Te, ne)};
		double ln_alpha {std::log(debye_length / expect_b0)};
		if (std::isnan(ln_alpha))
		{
			#pragma omp critical
			{
				std::cerr << "Error! ln_alpha = nan\n";
				std::cerr << "  Te = " << Te << '\n';
				std::cerr << "  ne = " << ne << '\n';
				std::cerr << "  charge = " << imp.get_charge() << '\n';
				std::cerr << "  debye_length = " << debye_length << '\n';
				std::cerr << "  expect_b0 = " << expect_b0 << '\n';
				std::cerr << "  ln_alpha = " << ln_alpha << '\n';
			}
		}

		// Calculate s (Eq. 19). Assume ni=ne and singly charge background
		double square_term {imp.get_charge() * Constants::charge_e 
			* Constants::charge_e / (Constants::eps0 * mu_ab)};
		double s {ln_alpha / (4.0 * Constants::pi) * square_term * square_term
			* ne / (g*g*g) * imp_time_step};

		// I believe this should be avoided now, but I leave this error message
		// since you never know.
		if (std::isnan(s))
		{
			#pragma omp critical
			{
				std::cerr << "Error! s = nan\n";
				std::cerr << "  T = " << T << '\n';
				std::cerr << "  Te = " << Te << '\n';
				std::cerr << "  bkg_vX = " << bkg_vX << '\n';
				std::cerr << "  bkg_vY = " << bkg_vY << '\n';
				std::cerr << "  bkg_vZ = " << bkg_vZ << '\n';
				std::cerr << "  imp_vX = " << imp.get_vX() << '\n';
				std::cerr << "  imp_vY = " << imp.get_vY() << '\n';
				std::cerr << "  imp_vZ = " << imp.get_vZ() << '\n';
				std::cerr << "  g = " << g << '\n';
				std::cerr << "  expect_g_sq = " << expect_g_sq << '\n';
				std::cerr << "  ne = " << ne << '\n';
				std::cerr << "  charge = " << imp.get_charge() << '\n';
				std::cerr << "  debye_length = " << debye_length << '\n';
				std::cerr << "  expect_b0 = " << expect_b0 << '\n';
				std::cerr << "  ln_alpha = " << ln_alpha << '\n';
			}
		}

		return std::make_tuple(s, gX, gY, gZ);
	}

	// Calculate A term in Nanbu collision model
	double nanbu_calc_A(const double s)
	{
		// We apply the approximation from Nanbu that generally covers outside
		// the range of applicability.
		if (s < 0.01)
		{
			return 1.0 / s;
		}
		else if (s > 5.0)
		{
			return 3.0 * std::exp(-s);
		}
		else
		{
			// Use interpolation function to find A value from precomputed 
			// arrays of A(s) values found in nanbu_s_a.h. These values were 
			// calculated using python/calc_nanbu.py.
			return Utilities::linear_interpolate(Nanbu::s, Nanbu::A, s);
		}

	}

	// Calculate chi (deflection angle) in Nanbu collision model
	double nanbu_calc_chi(const double s, const double A)
	{
		// Random number uniformly distributed between 0-1
		double U {Random::get(0.0, 1.0)};

		// Calculate and return chi. If A is too large we can get an exponential
		// overflow. So in that case we apply the approximation. 
		double cos_chi {};
		if (A > 50.0)
		{
			cos_chi = 1.0 + s * std::log(U);
		}
		else
		{
			// If A is too small, then sinh can return nan. Physically, this
			// means the plasma is too collisionless, so there is no real
			// deflection, so we return 0 if that happens. 
			if (A < 1e-5) return 0.0;

			// Normal calculation
			cos_chi = 1.0 / A * std::log(std::exp(-A) + 2.0 * U * std::sinh(A));
		}
		return std::acos(cos_chi);

	}

	// Calculate post-collision Cartesian velocity components from Nanbu model
	std::tuple<double, double, double> nanbu_post_coll(const double gX,
		const double gY, const double gZ, const double chi, const double mass_a,
		const double mass_b, const Impurity::Impurity& imp)
	{
		// For calculating h components
		double g_perp {std::sqrt(gY*gY + gZ*gZ)};
		double g {std::sqrt(gX*gX + gY*gY + gZ*gZ)};
		double eps {2.0 * Constants::pi * Random::get(0.0, 1.0)};
		double cos_eps {std::cos(eps)};
		double sin_eps {std::sin(eps)};

		// Calculate the h components
		double hX {g_perp * cos_eps};
		double hY {-(gY * gX * cos_eps + g * gZ * sin_eps) / g_perp};
		double hZ {-(gZ * gX * cos_eps - g * gY * sin_eps) / g_perp};

		// Calculate post collision velocities
		double mu {mass_b / (mass_a + mass_b)};
		double cos_chi {std::cos(chi)};
		double sin_chi {std::sin(chi)};
		double vX_post {imp.get_vX() - mu * (gX * (1.0 - cos_chi) 
			+ hX * sin_chi)};
		double vY_post {imp.get_vY() - mu * (gY * (1.0 - cos_chi) 
			+ hY * sin_chi)};
		double vZ_post {imp.get_vZ() - mu * (gZ * (1.0 - cos_chi) 
			+ hZ * sin_chi)};

		// Return as tuple
		return std::make_tuple(vX_post, vY_post, vZ_post);
	}

	void nanbu_coll(Impurity::Impurity& imp, const Background::Background& bkg,
		const int tidx, const int xidx, const int yidx, const int zidx,
		const Options::Options& opts, bool elec, const double imp_time_step)
	{
		// Derivation and steps taken from:
		// Nanbu, K. Theory of cumulative small-angle collisions in plasmas. 
		// Phys. Rev. E 55, 4642–4652 (1997).

		// A neutral will not experience a Coloumb collision (in fact, will
		// cause a divide by zero later on in this algorithm), so do nothing
		// in that case.
		if (imp.get_charge() == 0) return;

		// Will always need electron temperature/density
		double Te {bkg.get_te()(tidx, xidx, yidx, zidx)};
		double ne {bkg.get_ne()(tidx, xidx, yidx, zidx)};

		// Load some reusable variables based on which species. If elec = true,
		// then electrons and ions if not. 
		double mass_kg {};
		double T {};
		if (elec)
		{
			mass_kg = opts.gkyl_elec_mass_amu() * Constants::amu_to_kg;	
			T = Te;
		}
		else
		{
			mass_kg = opts.gkyl_ion_mass_amu() * Constants::amu_to_kg;	
			T = bkg.get_ti()(tidx, xidx, yidx, zidx);
		}

		// Calculate s (Eq. 19)
		auto [s, gX, gY, gZ] = nanbu_calc_s(imp, bkg, tidx, xidx, yidx, zidx, 
			imp_time_step, T, mass_kg, Te, ne);

		// Calcluate A (Eq. 13)
		double A {nanbu_calc_A(s)};

		// Calculate deflection angle, chi (Eq. 17)
		double chi {nanbu_calc_chi(s, A)};
		if (std::isnan(chi))
		{
			std::cerr << "Error! chi = nan\n";
			std::cerr << "  T = " << T << '\n';
			std::cerr << "  Te = " << Te << '\n';
			std::cerr << "  ne = " << ne << '\n';
			std::cerr << "  chi = " << chi << '\n';
			std::cerr << "  s = " << s << '\n';
			std::cerr << "  A = " << A << '\n';
		}

		// Calculate post-collision velocity (Eq. 20a)
		auto [vX_post, vY_post, vZ_post] = nanbu_post_coll(gX, gY, gZ, chi, 
			imp.get_mass(), mass_kg, imp);
		//std::cout << "pre: vX, vY, vZ = " << imp.get_vX() << ", " 
		//	<< imp.get_vY() << ", " << imp.get_vZ() << '\n';
		//std::cout << "post: vX, vY, vZ = " << vX_post << ", " << vY_post << ", " 
		//	<< vZ_post << '\n';

		// Update impurity
		imp.set_vX(vX_post);
		imp.set_vY(vY_post);
		imp.set_vZ(vZ_post);


	}
}
