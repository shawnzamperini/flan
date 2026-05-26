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
#include "impurity_stats.h"
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
		const Background::Background& bkg, Impurity::Impurity& imp, 
		const double T, const double uX, const double uY, const double uZ,
		const double m)
	{
		// Sampling from a drifting Maxwellian, which can be done by just
		// sampling from a normally distributed number with uXYZ = mean = drift 
		// and mu = sqrt(kT/m). We call this thermal sampling of the background
		// velocity. Unfortunately throwing away the second random number 
		// here. T [eV], m [kg]
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
		const double imp_time_step, const double T, const double mass_kg, 
		const double Te, const double ne, const double imp_vX_n, 
		const double imp_vY_n, const double imp_vZ_n)
	{
		// While I think the below is still valid, it was solved in calc_chi
		// by returning an isotropically distributed chi at large s (small A)
		// values. Can likely remove this part of things.
		// ---
		// There is an insanely subtle thing to comment on here that has real
		// effects. The Nanbu algorithm uses the instantaneous g to calcuate
		// s. This means that particle velocities far from the background
		// velocity get disaproportionally large kicks since s ~ g-3. This
		// can cause the average many-particle velocity to actually overshoot
		// the mean flow as they accelerate up to it from rest, ignoring all
		// other collisional phenomena. This is technically incorrect and 
		// unphysical, so we modify the algorithm as such:
		//   - Sample the background velocity including thermal sampling,
		//     vX, vY, vZ. This determines the X,Y,Z direction of g and allows
		//     random collisions with randomly directed background particles
		//	   to stay part of the algorithm. We call these the instantaneous
		//     relative velocities.
		//   - Enforce |g| to be that from using the flow velocity, no thermal
		//     sampling, and rescale gX, gY, gZ. This prevents overshooting 
		//     the mean flow. We call this the mean relative velocity.

		// Background flow at impurity location
		double uX {bkg.interp_uX_at_imp(imp)};
		double uY {bkg.interp_uY_at_imp(imp)};
		double uZ {bkg.interp_uZ_at_imp(imp)};

		// Random sample of background species instantanoues velocity 
		// (flow + thermal sampling) that particle is colliding with.
		auto [bkg_vX, bkg_vY, bkg_vZ] = sample_bkg_velocity(bkg, imp, T, uX, 
			uY, uZ, mass_kg);

		// XYZ components of instantaneous relative velocity
		double inst_gX {imp.get_vX() - bkg_vX};
		double inst_gY {imp.get_vY() - bkg_vY};
		double inst_gZ {imp.get_vZ() - bkg_vZ};

		// XYZ components of mean relative velocity
		double mean_gX {imp.get_vX() - uX};
		double mean_gY {imp.get_vY() - uY};
		double mean_gZ {imp.get_vZ() - uZ};

		// Magnitude of instantaneous and mean relative velocities.
		double mean_g {std::sqrt(mean_gX*mean_gX + mean_gY*mean_gY 
			+ mean_gZ*mean_gZ)};
		double inst_g {std::sqrt(inst_gX*inst_gX + inst_gY*inst_gY 
			+ inst_gZ*inst_gZ)};

		// Limit to a sufficiently small number to avoid overflows that can 
		// occur in std::exp and std::sinh later on.
		mean_g = std::max(mean_g, 0.001);
		inst_g = std::max(inst_g, 0.001);

		// Scale the instantanous velocity magnitude to equal the mean. 
		// Conceptually there is a basis for this, but I am not 100% convinced
		// and I'd rather just leave the algorithm alone and not modify it if
		// not necessary.
		//inst_gX = inst_gX * mean_g / inst_g;
		//inst_gY = inst_gY * mean_g / inst_g;
		//inst_gZ = inst_gZ * mean_g / inst_g;

		// Reduced mass of impurity and species colliding with
		double mu_ab {imp.get_mass() * mass_kg / (imp.get_mass() + mass_kg)};

		// Expectation value of g^2. Limit to a sufficiently small number to
		// avoid overflows that can occur in std::exp and std::sinh.
		double expect_g_sq {3.0 * T * Constants::ev_to_j / mu_ab};
		expect_g_sq = std::max(expect_g_sq, 0.001);

		// Expectation value of b0. Assuming singly charged background species.
		double expect_b0 {imp.get_charge() * Constants::charge_e 
			* Constants::charge_e / (2.0 * Constants::pi * Constants::eps0 
			* mu_ab * expect_g_sq)};

		// Debye length and Coulumb logarithm
		double debye_length {calc_debye_length(Te, ne)};
		double ln_alpha {std::log(debye_length / expect_b0)};

		// We can't let debye_length < b0, this can mess up the Nanbu model
		// because it makes ln_alpha < 0. In real plasmas it shouldn't go 
		// below 1, so clamp it.
		ln_alpha = std::max(ln_alpha, 1.0);

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

		// Calculate s (Eq. 19). Assume ni=ne and singly charge background.
		double square_term {imp.get_charge() * Constants::charge_e 
			* Constants::charge_e / (Constants::eps0 * mu_ab)};

		// The correct and original kinetic calculaton of s
		double s {ln_alpha / (4.0 * Constants::pi) * square_term * square_term
			* ne / (inst_g*inst_g*inst_g) * imp_time_step};

		// Needed to match assumptions made in the equation for the friction 
		// force
		//double s {ln_alpha / (4.0 * Constants::pi) * square_term * square_term
		//	* ne / (mean_g*mean_g*mean_g) * imp_time_step};
			
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
				std::cerr << "  imp_vX_n = " << imp_vX_n << '\n';
				std::cerr << "  imp_vY_n = " << imp_vY_n << '\n';
				std::cerr << "  imp_vZ_n = " << imp_vZ_n << '\n';
				std::cerr << "  imp_vX = " << imp.get_vX() << '\n';
				std::cerr << "  imp_vY = " << imp.get_vY() << '\n';
				std::cerr << "  imp_vZ = " << imp.get_vZ() << '\n';
				std::cerr << "  inst_g = " << inst_g << '\n';
				std::cerr << "  inst_gX = " << inst_gX << '\n';
				std::cerr << "  inst_gY = " << inst_gY << '\n';
				std::cerr << "  inst_gZ = " << inst_gZ << '\n';
				//std::cerr << "  mean_g = " << mean_g << '\n';
				//std::cerr << "  mean_gX = " << mean_gX << '\n';
				//std::cerr << "  mean_gY = " << mean_gY << '\n';
				//std::cerr << "  mean_gZ = " << mean_gZ << '\n';
				std::cerr << "  expect_g_sq = " << expect_g_sq << '\n';
				std::cerr << "  ne = " << ne << '\n';
				std::cerr << "  charge = " << imp.get_charge() << '\n';
				std::cerr << "  debye_length = " << debye_length << '\n';
				std::cerr << "  expect_b0 = " << expect_b0 << '\n';
				std::cerr << "  ln_alpha = " << ln_alpha << '\n';
			}
		}

		return std::make_tuple(s, inst_gX, inst_gY, inst_gZ);
	}

	// Calculate A term in Nanbu collision model
	double nanbu_calc_A(const double s)
	{
		// Not very collisional approximation
		if (s < 0.01)
		{
			return 1.0 / s;
		}

		// Very collisional approximation
		else if (s > 5.0)
		{
			return 3.0 * std::exp(-s);
		}

		// Little bit of collisional, little bit of not collisional
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

		// Boolean masks to determine which case we want (small, middle or 
		// large s calculation). This avoids branching and makes for for
		// SIMD-friendly code for the compiler to vectorize.
		const bool mask_small {(s < 0.02)};
		const bool mask_large {(s > 6.0)};

		// Compute cos_chi when s is small or large
		const double cos_small {1.0 + s * std::log(U)};
		const double cos_large {2.0 * U - 1.0};

		// Compute cos_chi when s is somewhere in between
		const double expA  {std::exp(-A)};
		const double sinhA {std::sinh(A)};
		const double cos_mid {(1.0 / A) * std::log(expA + 2.0 * U * sinhA)};

		// Blend to choose the correct branch
		double cos_chi {cos_mid};
		cos_chi = mask_small ? cos_small : cos_chi;
		cos_chi = mask_large ? cos_large : cos_chi;

		return std::acos(cos_chi);
	}

	// Calculate post-collision Cartesian velocity components from Nanbu model
	std::tuple<double, double, double> nanbu_post_coll(const double gX,
		const double gY, const double gZ, const double chi, const double mass_a,
		const double mass_b, const Impurity::Impurity& imp, 
		const double imp_vX_n, const double imp_vY_n, const double imp_vZ_n)
	{
		// For calculating h components
		double g_perp {std::sqrt(gY*gY + gZ*gZ)};
		double g {std::sqrt(gX*gX + gY*gY + gZ*gZ)};
		double eps {2.0 * Constants::pi * Random::get(0.0, 1.0)};
		double cos_eps {std::cos(eps)};
		double sin_eps {std::sin(eps)};

		// Calculate the h components. If g_perp = 0 then g is
		// parallel to X and can cause NaNs/infs when we divide by g_perp, so
		// we handle accordingly.
		double hX {};
		double hY {};
		double hZ {};
		if (g_perp < Constants::small)
		{
			hX = 0.0;
			hY = g * cos_eps;
			hZ = g * sin_eps;
		}
		else
		{
			hX = g_perp * cos_eps;
			hY = -(gY * gX * cos_eps + g * gZ * sin_eps) / g_perp;
			hZ = -(gZ * gX * cos_eps - g * gY * sin_eps) / g_perp;
		}

		// Calculate post collision velocities
		double mu {mass_b / (mass_a + mass_b)};
		double cos_chi {std::cos(chi)};
		double sin_chi {std::sin(chi)};
		
		// From Nanbu paper, assuming self-consistent particle interaction
		double vX_post {imp.get_vX() - mu * (gX * (1.0 - cos_chi) 
			+ hX * sin_chi)};
		double vY_post {imp.get_vY() - mu * (gY * (1.0 - cos_chi) 
			+ hY * sin_chi)};
		double vZ_post {imp.get_vZ() - mu * (gZ * (1.0 - cos_chi) 
			+ hZ * sin_chi)};

		// AI suggested the sign in front of the h term was wrong
		//double vX_post {imp.get_vX() - mu * (gX * (1.0 - cos_chi) 
		//	- hX * sin_chi)};
		//double vY_post {imp.get_vY() - mu * (gY * (1.0 - cos_chi) 
		//	- hY * sin_chi)};
		//double vZ_post {imp.get_vZ() - mu * (gZ * (1.0 - cos_chi) 
		//	- hZ * sin_chi)};

		// From AI. It suggested that a) the sign infront of the h term was
		// wrong and b) that since we are in the test particle limit we don't
		// want mu and gave this derivation instead.
		//double vX_post {imp.get_vX() + gX * (cos_chi - 1.0) + hX * sin_chi};
		//double vY_post {imp.get_vY() + gY * (cos_chi - 1.0) + hY * sin_chi};
		//double vZ_post {imp.get_vZ() + gZ * (cos_chi - 1.0) + hZ * sin_chi};

		// Full time step option, unsure if needed...
		//double vX_post {imp_vX_n - mu * (gX * (1.0 - cos_chi) 
		//	+ hX * sin_chi)};
		//double vY_post {imp_vY_n - mu * (gY * (1.0 - cos_chi) 
		//	+ hY * sin_chi)};
		//double vZ_post {imp_vZ_n - mu * (gZ * (1.0 - cos_chi) 
		//	+ hZ * sin_chi)};

		// Return as tuple
		return std::make_tuple(vX_post, vY_post, vZ_post);
	}

	void nanbu_coll(Impurity::Impurity& imp, const Background::Background& bkg,
		const int tidx, const int xidx, const int yidx, const int zidx,
		const Options::Options& opts, bool elec, const double imp_time_step,
		Impurity::Statistics& imp_stats)
	{
		// Derivation and steps taken from:
		// Nanbu, K. Theory of cumulative small-angle collisions in plasmas. 
		// Phys. Rev. E 55, 4642–4652 (1997).

		// A neutral will not experience a Coloumb collision (in fact, will
		// cause a divide by zero later on in this algorithm), so do nothing
		// in that case.
		if (imp.get_charge() == 0) return;

		// Will always need electron temperature/density. Trilinearly
		// interpolate in space and then linearly interpolate in time.
		//double Te {bkg.get_te()(tidx, xidx, yidx, zidx)};
		//double ne {bkg.get_ne()(tidx, xidx, yidx, zidx)};
		double Te {bkg.interp_te_at_imp(imp)};
		double ne {bkg.interp_ne_at_imp(imp)};

		// Need to account for this better, just putting it here for now so I
		// can get this paper submitted :(
		if (ne < 1e16) ne = 1e16;
		if (Te < 0.1) ne = 0.1;

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
			//T = bkg.get_ti()(tidx, xidx, yidx, zidx);
			T = bkg.interp_ti_at_imp(imp);
		}

		// Subtlety! Impurity velocity is defined at half time steps, but the collision
		// happens at a full time step. Since this collision update is happening
		// after the Boris update, imp.v is at n+1/2 and imp.prev_v is at n-1/2.
		// So reconstruct the full-time step as
		// the average between v_n-1/2 and v_n+1/2 and use that in the following.
		// We will overwrite v_n+1/2 with the post-collision value.
		// ---
		// These aren't used right now in the following functions, but it still 
		// may be the correct thing to do. Will need time to test/compare the 
		// two. Right now just set to zero so we can leave the function calls
		// in place.
		//double imp_vX_n {(imp.get_vX() + imp.get_prev_vX()) / 2.0};
		//double imp_vY_n {(imp.get_vY() + imp.get_prev_vY()) / 2.0};
		//double imp_vZ_n {(imp.get_vZ() + imp.get_prev_vZ()) / 2.0};
		constexpr double imp_vX_n {0.0};
		constexpr double imp_vY_n {0.0};
		constexpr double imp_vZ_n {0.0};

		// The Nanbu model has three main variables in it:
		// s:    How collisional is this step?
		// A(s): What is the shape of the scattering distribution? A(s) is
		//       the PDF.
		// chi:  What is the actual deflection angle this time? It is calculated
		//       via direct inversion from the CDF formed from our PDF (A(s)). 

		// Calculate s (Eq. 19), making sure to pass in the full time step
		// velocities.
		auto [s, gX, gY, gZ] = nanbu_calc_s(imp, bkg, imp_time_step, T, 
			mass_kg, Te, ne, imp_vX_n, imp_vY_n, imp_vZ_n);

		// Add to running sum of s values in each cell so we can do an average
		// later. Only consider for ions since they are the dominant collision
		// but no reason this can't be expanded for electrons as well. 
		// Generally leave this commented out unless you are investigating the 
		// collision model.
		if (!elec)
		{
			imp_stats.add_s(tidx, xidx, yidx, zidx, static_cast<BkgFPType>(s));
		}

		// Calcluate A (Eq. 13)
		double A {nanbu_calc_A(s)};
		//std::cout << "s = " << s << "\tA = " << A << '\n';

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
			imp.get_mass(), mass_kg, imp, imp_vX_n, imp_vY_n, imp_vZ_n);
		//std::cout << "pre: vX, vY, vZ = " << imp.get_vX() << ", " 
		//	<< imp.get_vY() << ", " << imp.get_vZ() << '\n';
		//std::cout << "post: vX, vY, vZ = " << vX_post << ", " << vY_post << ", " 
		//	<< vZ_post << '\n';

		// Update impurity velocity. We have calculated the post-collision
		// v at the full time step n, v_n, but imp.v is v_n+1/2 due to the way
		// the Boris algorithm uses a leap-frog scheme. So we need to
		// accelerate the velocity forward half a timestep before storing it.
		// The Lorentz force is used to calculate the acceleration at time
		// step n.
		//const double q_m {imp.get_charge() * Constants::charge_e 
		//	/ imp.get_mass()};
		//const double EX {bkg.get_eX()(tidx, xidx, yidx, zidx)};
		//const double EY {bkg.get_eY()(tidx, xidx, yidx, zidx)};
		//const double EZ {bkg.get_eZ()(tidx, xidx, yidx, zidx)};
		//const double BX {bkg.get_bX()(tidx, xidx, yidx, zidx)};
		//const double BY {bkg.get_bY()(tidx, xidx, yidx, zidx)};
		//const double BZ {bkg.get_bZ()(tidx, xidx, yidx, zidx)};

		// Assemble arrays and calculate v x B
		//const std::array<double, 3> v_n {vX_post, vY_post, vZ_post};
		//const std::array<double, 3> B_n {BX, BY, BZ};
		//const std::array<double, 3> v_cross_B	
		//	{Utilities::cross_product(v_n, B_n)};

		// Calculate acceleration. This is just a = F_Lorentz / m
		//double aX_n {q_m * (EX + v_cross_B[0])};
		//double aY_n {q_m * (EY + v_cross_B[1])};
		//double aZ_n {q_m * (EZ + v_cross_B[2])};

		//imp.set_vX(vX_post + 0.5 * aX_n * imp_time_step, false);
		//imp.set_vY(vY_post + 0.5 * aY_n * imp_time_step, false);
		//imp.set_vZ(vZ_post + 0.5 * aZ_n * imp_time_step, false);

		// Update impurity velocity, which is v_n+1/2 since this is happening
		// after the Boris update.
		imp.set_vX(vX_post);
		imp.set_vY(vY_post);
		imp.set_vZ(vZ_post);

	}
}
