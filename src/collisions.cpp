/**
* @file  collisions.cpp
*
* @brief Routines for updating particle velocity due to Coulomb collisions via
*        the Nanbu algorithm.
*
* Derivation and steps taken from:
* Nanbu, K. Theory of cumulative small-angle collisions in plasmas.
* Phys. Rev. E 55, 4642–4652 (1997).
*
* Each time step, the cumulative effect of all the small-angle collisions
* with a background species is replaced by a single deflection of the
* relative velocity by an angle chi. How collisional the step is is set by
* s, which sets the shape A(s) of the distribution chi is sampled from.
*/
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <utility>
#include <vector>

#include "background.h"
#include "constants.h"
#include "impurity_stats.h"
#include "nanbu_s_a.h"
#include "options.h"
#include "pcg32.h"
#include "slots.h"
#include "utilities.h"
#include "variance_reduction.h"

namespace Collisions
{
	/**
	* Calculates the electron Debye length.
	*
	* @param te Electron temperature (eV).
	* @param ne Electron density (m^-3).
	* @return Debye length (m).
	*/
	double calc_debye_length(const double te, const double ne)
	{
		return std::sqrt(Constants::eps0 * te * Constants::ev_to_j / 
			(ne * Constants::charge_e * Constants::charge_e)); 
	}


	/**
	* Samples the velocity of a background particle from a drifting
	* Maxwellian.
	*
	* Each component is drawn from a normal distribution with a mean equal
	* to the flow velocity and a standard deviation of sqrt(kT/m). We call
	* this thermal sampling of the background velocity.
	*
	* @param T Temperature of the background species (eV).
	* @param uX Cartesian X flow velocity of the background species (m/s).
	* @param uY Cartesian Y flow velocity of the background species (m/s).
	* @param uZ Cartesian Z flow velocity of the background species (m/s).
	* @param m Mass of the background species (kg).
	* @param rng Random number generator for this thread.
	* @return Sampled velocity {vX, vY, vZ} in Cartesian space (m/s).
	*/
	std::tuple<double, double, double> sample_bkg_velocity(
		const double T, const double uX, const double uY, const double uZ,
		const double m, pcg32& rng)
	{
		// Sampling from a drifting Maxwellian, which can be done by just
		// sampling from a normally distributed number with uXYZ = mean = drift 
		// and mu = sqrt(kT/m). We call this thermal sampling of the background
		// velocity. Unfortunately throwing away the second random number 
		// here. T [eV], m [kg]
		double mu {std::sqrt(T * Constants::ev_to_j / m)};  // m/s
		double vX {rng.normal(uX, mu)};
		double vY {rng.normal(uY, mu)};
		double vZ {rng.normal(uZ, mu)};

		// Return as tuple
		return std::make_tuple(vX, vY, vZ);
	}


	/**
	* Calculates the s term in the Nanbu collision model (Eq. 19), which
	* measures how collisional this time step is.
	*
	* The relative velocity g is between the impurity and a background
	* particle whose velocity is sampled from a drifting Maxwellian (flow
	* plus thermal sampling) with sample_bkg_velocity. The Coulomb logarithm
	* uses the Debye length and the expectation value of the impact
	* parameter b0, and is clamped to be at least 1. The background species
	* is assumed to be singly charged with ni = ne. An error is printed if
	* the Coulomb logarithm or s is NaN.
	*
	* @param bkg Background plasma, used to interpolate the flow velocity.
	* @param imp_time_step Time step (s).
	* @param T Temperature of the background species colliding with the
	*        impurity (eV).
	* @param mass_kg Mass of the background species (kg).
	* @param Te Electron temperature (eV), used for the Debye length.
	* @param ne Electron density (m^-3).
	* @param t Time of the particle (s).
	* @param x Computational x coordinate of the particle.
	* @param y Computational y coordinate of the particle.
	* @param z Computational z coordinate of the particle.
	* @param vX Cartesian X velocity of the impurity (m/s).
	* @param vY Cartesian Y velocity of the impurity (m/s).
	* @param vZ Cartesian Z velocity of the impurity (m/s).
	* @param q Charge state of the impurity. Must not be 0.
	* @param slots Particle slots, used for the impurity mass.
	* @param rng Random number generator for this thread.
	* @return {s, gX, gY, gZ}, where gX, gY and gZ are the Cartesian
	*         components of the relative velocity (m/s).
	*/
	std::tuple<double, double, double, double> nanbu_calc_s(
		const Background::Background& bkg, 
		const double imp_time_step, const double T, const double mass_kg, 
		const double Te, const double ne, const double t,
		const double x, const double y, const double z, const double vX,
		const double vY, const double vZ, const int q, Slots::Slots& slots,
		pcg32& rng)
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
		double uX {bkg.interp_uX(t, x, y, z)}; 
		double uY {bkg.interp_uY(t, x, y, z)}; 
		double uZ {bkg.interp_uZ(t, x, y, z)}; 

		// Random sample of background species instantanoues velocity 
		// (flow + thermal sampling) that particle is colliding with.
		auto [bkg_vX, bkg_vY, bkg_vZ] = sample_bkg_velocity(T, uX, 
			uY, uZ, mass_kg, rng);

		// XYZ components of instantaneous relative velocity
		double inst_gX {vX - bkg_vX};
		double inst_gY {vY - bkg_vY};
		double inst_gZ {vZ - bkg_vZ};

		// XYZ components of mean relative velocity
		//double mean_gX {vX - uX};
		//double mean_gY {vY - uY};
		//double mean_gZ {vZ - uZ};

		// Magnitude of instantaneous and mean relative velocities.
		//double mean_g {std::sqrt(mean_gX*mean_gX + mean_gY*mean_gY 
		//	+ mean_gZ*mean_gZ)};
		double inst_g {std::sqrt(inst_gX*inst_gX + inst_gY*inst_gY 
			+ inst_gZ*inst_gZ)};

		// Limit to a sufficiently small number to avoid overflows that can 
		// occur in std::exp and std::sinh later on.
		//mean_g = std::max(mean_g, 0.001);
		inst_g = std::max(inst_g, 0.001);

		// Scale the instantanous velocity magnitude to equal the mean. 
		// Conceptually there is a basis for this, but I am not 100% convinced
		// and I'd rather just leave the algorithm alone and not modify it if
		// not necessary.
		//inst_gX = inst_gX * mean_g / inst_g;
		//inst_gY = inst_gY * mean_g / inst_g;
		//inst_gZ = inst_gZ * mean_g / inst_g;

		// Reduced mass of impurity and species colliding with
		double mu_ab {slots.mass() * mass_kg / (slots.mass() + mass_kg)};

		// Expectation value of g^2. Limit to a sufficiently small number to
		// avoid overflows that can occur in std::exp and std::sinh.
		double expect_g_sq {3.0 * T * Constants::ev_to_j / mu_ab};
		expect_g_sq = std::max(expect_g_sq, 0.001);

		// Expectation value of b0. Assuming singly charged background species.
		double expect_b0 {q * Constants::charge_e 
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
				std::cerr << "  charge = " << q << '\n';
				std::cerr << "  debye_length = " << debye_length << '\n';
				std::cerr << "  expect_b0 = " << expect_b0 << '\n';
				std::cerr << "  ln_alpha = " << ln_alpha << '\n';
			}
		}

		// Calculate s (Eq. 19). Assume ni=ne and singly charge background.
		double square_term {q * Constants::charge_e 
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
				std::cerr << "  s = " << s << '\n';
				std::cerr << "  T = " << T << '\n';
				std::cerr << "  Te = " << Te << '\n';
				std::cerr << "  bkg_vX = " << bkg_vX << '\n';
				std::cerr << "  bkg_vY = " << bkg_vY << '\n';
				std::cerr << "  bkg_vZ = " << bkg_vZ << '\n';
				std::cerr << "  imp_vX = " << vX << '\n';
				std::cerr << "  imp_vY = " << vY << '\n';
				std::cerr << "  imp_vZ = " << vZ << '\n';
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
				std::cerr << "  charge = " << q << '\n';
				std::cerr << "  debye_length = " << debye_length << '\n';
				std::cerr << "  expect_b0 = " << expect_b0 << '\n';
				std::cerr << "  ln_alpha = " << ln_alpha << '\n';
			}
		}

		return std::make_tuple(s, inst_gX, inst_gY, inst_gZ);
	}


	/**
	* Calculates the A term in the Nanbu collision model (Eq. 13), which sets
	* the shape of the distribution the deflection angle is sampled from.
	*
	* For s < 0.01 (weakly collisional) A = 1 / s, and for s > 5 (strongly
	* collisional) A = 3 exp(-s). In between, A is linearly interpolated
	* from the precomputed table in nanbu_s_a.h, which was calculated using
	* python/calc_nanbu.py.
	*
	* @param s The s term from nanbu_calc_s.
	* @return The A term.
	*/
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
			// We use the uniform interpolate since it's a lot faster, just
			// kinda leaving this comment here to point that out.
			//return Utilities::linear_interpolate(Nanbu::s, Nanbu::A, s);
			return Utilities::linear_interpolate_uniform(Nanbu::s, Nanbu::A, s);
		}

	}


	/**
	* Samples the deflection angle chi in the Nanbu collision model (Eq. 17).
	*
	* cos(chi) is found by directly inverting the cumulative distribution
	* built from A, using a uniform random number U. For s < 0.02 the small-s
	* limit 1 + s ln(U) is used, and for s > 6 the deflection is isotropic,
	* 2U - 1. cos(chi) is clamped to [-1, 1] to guard against rounding
	* errors.
	*
	* @param s The s term from nanbu_calc_s.
	* @param A The A term from nanbu_calc_A.
	* @param rng Random number generator for this thread.
	* @return Deflection angle chi (radians).
	*/
	double nanbu_calc_chi(const double s, const double A, pcg32& rng)
	{
		// Random number uniformly distributed between 0-1
		double U {rng.next_double()};
		U = std::max(U, 1e-12);

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

		// Clamp to handle any rounding errors that could error acos
		cos_chi = std::clamp(cos_chi, -1.0, 1.0);
		return std::acos(cos_chi);
	}


	/**
	* Calculates the post-collision Cartesian velocity of the impurity in the
	* Nanbu collision model (Eq. 20a).
	*
	* The relative velocity g is rotated by the deflection angle chi about a
	* random azimuthal angle, uniformly distributed in [0, 2 pi). The case
	* where g is parallel to X is handled separately to avoid dividing by
	* zero.
	*
	* @param gX Cartesian X component of the relative velocity (m/s).
	* @param gY Cartesian Y component of the relative velocity (m/s).
	* @param gZ Cartesian Z component of the relative velocity (m/s).
	* @param chi Deflection angle from nanbu_calc_chi (radians).
	* @param mass_a Mass of the impurity (kg).
	* @param mass_b Mass of the background species (kg).
	* @param vX Cartesian X velocity of the impurity before the collision
	*        (m/s).
	* @param vY Cartesian Y velocity of the impurity before the collision
	*        (m/s).
	* @param vZ Cartesian Z velocity of the impurity before the collision
	*        (m/s).
	* @param rng Random number generator for this thread.
	* @return Post-collision velocity {vX, vY, vZ} in Cartesian space (m/s).
	*/
	std::tuple<double, double, double> nanbu_post_coll(const double gX,
		const double gY, const double gZ, const double chi, const double mass_a,
		const double mass_b, const double vX, const double vY, const double vZ,
		pcg32& rng)
	{
		// For calculating h components
		double g_perp {std::sqrt(gY*gY + gZ*gZ)};
		double g {std::sqrt(gX*gX + gY*gY + gZ*gZ)};
		double eps {2.0 * Constants::pi * rng.next_double()};
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
		double vX_post {vX - mu * (gX * (1.0 - cos_chi) + hX * sin_chi)};
		double vY_post {vY - mu * (gY * (1.0 - cos_chi) + hY * sin_chi)};
		double vZ_post {vZ - mu * (gZ * (1.0 - cos_chi) + hZ * sin_chi)};

		// Return as tuple
		return std::make_tuple(vX_post, vY_post, vZ_post);
	}


	/**
	* Updates the Cartesian velocity of every particle due to Coulomb
	* collisions with one background species using the Nanbu model.
	*
	* Dead and neutral particles are skipped. For each particle, ne and Te
	* (and Ti for ions) are interpolated at the particle location, then s,
	* A and chi are calculated and the velocity is updated with
	* nanbu_post_coll. Since this happens after the Boris update, the
	* updated velocity is v at n+1/2. An error is printed if ne < 0 or chi is
	* NaN.
	*
	* For ion collisions, s * weight * dt is also added to imp_stats so the
	* average s in each cell can be calculated later.
	*
	* @param slots Particle slots. The velocities are updated in place.
	* @param bkg Background plasma.
	* @param opts Simulation options, used for the background masses
	*        (gkyl_elec_mass_amu, gkyl_ion_mass_amu) and min_ti.
	* @param elec If true, collide with electrons (Te). Otherwise collide
	*        with ions (Ti, limited to at least min_ti).
	* @param dt Time step (s).
	* @param imp_stats Impurity statistics, where the s values are summed
	*        for ion collisions.
	* @param rngs Random number generators, one per OpenMP thread.
	*/
	void nanbu_coll(Slots::Slots& slots, const Background::Background& bkg,
		const Options::Options& opts, bool elec, const double dt, 
		Impurity::Statistics& imp_stats, std::vector<pcg32>& rngs)
	{
		// Derivation and steps taken from:
		// Nanbu, K. Theory of cumulative small-angle collisions in plasmas. 
		// Phys. Rev. E 55, 4642–4652 (1997).

		#pragma omp parallel
		{

			// Grab our RNG for this thread, each thread has it own (and they're
			// seeded uniquely).
			int tid = omp_get_thread_num();
			pcg32& rng = rngs[tid];

		#pragma omp for
		for (int i=0; i < slots.N(); ++i)
		{
			// Skip dead particles
			if (slots.state()[i] > 0) continue;
	
			// A neutral will not experience a Coloumb collision (in fact, will
			// cause a divide by zero later on in this algorithm), so do nothing
			// in that case.
			if (slots.q()[i] == 0) continue;

			// Local variables
			double t {slots.t()[i]};
			double x {slots.x()[i]};
			double y {slots.y()[i]};
			double z {slots.z()[i]};
			double vX {slots.vX()[i]};
			double vY {slots.vY()[i]};
			double vZ {slots.vZ()[i]};
			int state {slots.state()[i]};
			int q {slots.q()[i]};

			// Will always need electron temperature/density. Trilinearly
			// interpolate in space and then linearly interpolate in time.
			double ne {bkg.interp_ne(t, x, y, z)}; 
			double Te {bkg.interp_te(t, x, y, z)}; 

			// Don't let it go below a minimum. This can prevent extrapolating
			// near the edges going below zero.
			//if (ne < opts.min_ne())
			//{
			//	ne = std::max(ne, opts.min_ne());
			//}
			//Te = std::max(Te, opts.min_te());

	
			// This will ultimately cause a crash. I try to put in anough gaurds
			// to prevent this from happening (e.g., don't interpolate at the
			// edges of the grid), but it may be able to still slip in and
			// can be tricky to debug.
			if (ne < 0) 
			{
				#pragma omp critical
				{
					std::cout << "error! ne = " << ne << "\n";
					std::cout << "  txyz = " << t << ", " << x << ", " << y 
						<< ", " << z << "  state = " << state << '\n';
				}
			}

			// Don't let values go too low. Could happen with interpolation
			// near the edge of a grid, though we try to prevent that from
			// happening in the inteprolation functions.
			Te = std::max(opts.min_te(), Te);
			ne = std::max(opts.min_ne(), ne);

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
				T = bkg.interp_ti(t, x, y, z); 
				T = std::max(opts.min_ti(), T);
			}

			// The Nanbu model has three main variables in it:
			// s:    How collisional is this step?
			// A(s): What is the shape of the scattering distribution? A(s) is
			//       the PDF.
			// chi:  What is the actual deflection angle this time? It is calculated
			//       via direct inversion from the CDF formed from our PDF (A(s)). 

			// Calculate s (Eq. 19), making sure to pass in the full time step
			// velocities.
			auto [s, gX, gY, gZ] = nanbu_calc_s(bkg, dt, T, 
				mass_kg, Te, ne, t, x, y, z, vX, vY, vZ, q, slots, rng);

			// Add to running sum of s values in each cell so we can do an average
			// later. Only consider for ions since they are the dominant collision
			// but no reason this can't be expanded for electrons as well. 
			// Generally leave this commented out unless you are investigating the 
			// collision model.
			/*
			if (!elec)
			{
				double p_w {slots.weight()[i]};
				int tidx {slots.tidx()[i]};
				int xidx {slots.xidx()[i]};
				int yidx {slots.yidx()[i]};
				int zidx {slots.zidx()[i]};

				#pragma omp critical
				imp_stats.add_s(tidx, xidx, yidx, zidx, s * p_w * dt);
			}
			*/

			// Calcluate A (Eq. 13)
			double A {nanbu_calc_A(s)};

			// Calculate deflection angle, chi (Eq. 17)
			double chi {nanbu_calc_chi(s, A, rng)};
			if (std::isnan(chi))
			{
				#pragma omp critical
				{
					std::cerr << "Error! chi = nan\n";
					std::cerr << "  T = " << T << '\n';
					std::cerr << "  Te = " << Te << '\n';
					std::cerr << "  ne = " << ne << '\n';
					std::cerr << "  chi = " << chi << '\n';
					std::cerr << "  s = " << s << '\n';
					std::cerr << "  A = " << A << '\n';
				}
			}

			// Calculate post-collision velocity (Eq. 20a)
			auto [vX_post, vY_post, vZ_post] = nanbu_post_coll(gX, gY, gZ, chi, 
				slots.mass(), mass_kg, vX, vY, vZ, rng);

			// Update impurity velocity, which is v_n+1/2 since this is happening
			// after the Boris update.
			slots.set_vX(i, vX_post);
			slots.set_vY(i, vY_post);
			slots.set_vZ(i, vZ_post);

		}  // slots loop
		} // omp parallel
	}  // nanbu_coll

}  //namespace Collisions
