#include <cuda_runtime.h>
#include <cstdio>
#include <math.h>

#include "background_device.h"
#include "constants.h"
#include "pcg32.h"
#include "slots_device.h"


namespace Collisions
{
	__device__
	double calc_debye_length_cuda(const double te, const double ne)
	{
		return sqrt(Constants::eps0 * te * Constants::ev_to_j / 
			(ne * Constants::charge_e * Constants::charge_e)); 
	}


	__device__
	void sample_bkg_velocity(const double T, const double uX, const double uY, 
		const double uZ, const double m, pcg32& rng, double& bkg_vX, 
		double& bkg_vY, double& bkg_vZ)
	{
		// Sampling from a drifting Maxwellian, which can be done by just
		// sampling from a normally distributed number with uXYZ = mean = drift 
		// and mu = sqrt(kT/m). We call this thermal sampling of the background
		// velocity. Unfortunately throwing away the second random number 
		// here. T [eV], m [kg]
		double mu {std::sqrt(T * Constants::ev_to_j / m)};  // m/s
		bkg_vX = rng.normal(uX, mu);
		bkg_vY = rng.normal(uY, mu);
		bkg_vZ = rng.normal(uZ, mu);
	}


	// Calculate s term in Nanbu collision model
	__device__
	void nanbu_calc_s(const Background::BackgroundDevice& bkg_d, 
		const double imp_time_step, const double T, const double mass_kg, 
		const double Te, const double ne, const double t,
		const double x, const double y, const double z, const double vX,
		const double vY, const double vZ, const int q, 
		Slots::SlotsDevice& slots_d, pcg32& rng)
	{
		/*

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
		double mean_gX {vX - uX};
		double mean_gY {vY - uY};
		double mean_gZ {vZ - uZ};

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
				std::cerr << "  mean_g = " << mean_g << '\n';
				std::cerr << "  mean_gX = " << mean_gX << '\n';
				std::cerr << "  mean_gY = " << mean_gY << '\n';
				std::cerr << "  mean_gZ = " << mean_gZ << '\n';
				std::cerr << "  expect_g_sq = " << expect_g_sq << '\n';
				std::cerr << "  ne = " << ne << '\n';
				std::cerr << "  charge = " << q << '\n';
				std::cerr << "  debye_length = " << debye_length << '\n';
				std::cerr << "  expect_b0 = " << expect_b0 << '\n';
				std::cerr << "  ln_alpha = " << ln_alpha << '\n';
			}
		}

		return std::make_tuple(s, inst_gX, inst_gY, inst_gZ);
		*/
	}
}
