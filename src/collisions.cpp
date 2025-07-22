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
	double calc_reduced_mass(double m1, double m2)
	{
		return m1 * m2 / (m1 + m2);
	}

	double calc_debye_length(const double te, const double ne)
	{
		return std::sqrt(Constants::eps0 * te * Constants::ev_to_j / 
			(ne * Constants::charge_e * Constants::charge_e)); 
	}

	std::vector<double> generate_random_vels(double mu_vel, double sigma_vel, 
		int num)
	{
		// Loop through and fill out vector one random velocity at a time. We
		// use a while loop since we get two random numbers at a time.
		std::vector<double> vels (num);
		int count {};
		while (count < num)
		{
			// Get two normally distributed velocities and store into vels.
			std::pair<double, double> rans {Random::get_two_norm(mu_vel, 
				sigma_vel)};
			
			// If this returns a negative value just cap it at 0.0.
			if (rans.first < 0.0) vels[count] = 0.0;
			else vels[count] = rans.first;
			
			// Check that count + 1 does not go out of bounds. This can happen
			// if num is odd.
			if ((count + 1) < num){
				if (rans.second < 0.0) vels[count + 1] = 0.0;
				else vels[count + 1] = rans.second;
			}
			count += 2;
		}

		return vels;
	}

	// Calculate impurity velocity magnitude
	double calc_imp_vel(const Impurity::Impurity& imp)
	{
		return std::sqrt(imp.get_vX() * imp.get_vX() 
			+ imp.get_vY() * imp.get_vY() + imp.get_vZ() * imp.get_vZ());
	}
	
	std::array<double, 3> ran_unit_vector()
	{

		// Create empty vector, fill with random values between (-1,1) and
		// normalize. 
		std::array<double, 3> unit_vec {0.0, 0.0, 0.0};
		for (std::size_t i {}; i < 3; ++i) 
			unit_vec[i] = Random::get(-1.0, 1.0);
		double mag {};
		for (std::size_t i {}; i < 3; ++i) 
			mag += unit_vec[i] * unit_vec[i];
		mag = std::sqrt(mag);
		for (std::size_t i {}; i < 3; ++i) 
			unit_vec[i] /= mag;

		return unit_vec;
	}
	
	// Function to create a rotation matrix given an axis and angle. Thanks
	// to Copilot AI for this one!
	std::array<std::array<double, 3>, 3> calc_rotation_matrix(
		const std::array<double, 3>& axis, double angle) 
	{
		double x {axis[0]};
		double y {axis[1]};
		double z {axis[2]};
		double cosA {std::cos(angle)};
		double sinA {std::sin(angle)};
		double oneMinusCosA {1.0 - cosA};

		return {{
			{cosA + x*x*oneMinusCosA, x*y*oneMinusCosA - z*sinA, 
				x*z*oneMinusCosA + y*sinA},
			{y*x*oneMinusCosA + z*sinA, cosA + y*y*oneMinusCosA, 
				y*z*oneMinusCosA - x*sinA},
			{z*x*oneMinusCosA - y*sinA, z*y*oneMinusCosA + x*sinA, 
				cosA + z*z*oneMinusCosA}
		}};
	}
	
	// Function to multiply a vector by a matrix. Thanks to Copilot AI for 
	// this one!
	std::array<double, 3> multiply_matrix_vector(
		const std::array<std::array<double, 3>, 3>& matrix, 
		const std::array<double, 3>& vec) 
	{
		std::array<double, 3> result {0.0, 0.0, 0.0};
		for (int i = 0; i < 3; ++i) 
		{
			for (int j = 0; j < 3; ++j) 
			{
				result[i] += matrix[i][j] * vec[j];
			}
		}
		return result;
	}

	void rotate_imp(Impurity::Impurity& imp, double defl_ang)
	{
		// I had some help from Wikipedia and AI to figure this one out. We 
		// wish to create a new vector that is rotated by defl_ang in a random
		// direction. The steps to do this are:
		//  1. Create a random unit vector that we will rotate around.
		//  2. Compute the rotation matrix to rotate around the random unit 
		//     vector. See this section on Wikipedia for the matrix form:
		//     https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix
		//     _from_axis_and_angle
		//  3. Apply rotation by multiplying the original vector by the
		//     rotation matrix. 
		//  4. Assign new velocity components to impurity.

		// 1. Create random unit vector
		std::array<double, 3> rotation_axis {ran_unit_vector()};

		// 2. Calculate the rotation matrix. Don't forget to convert angle
		// to radians first.
		std::array<std::array<double, 3>, 3> rotation_matrix {
			calc_rotation_matrix(rotation_axis, 
			defl_ang)};

		// 3. Apply rotation matrix to find new rotated vector
		std::array<double, 3> init_vec {imp.get_vX(), imp.get_vY(), 
			imp.get_vZ()};
		std::array<double, 3> rotated_vec {multiply_matrix_vector(
			rotation_matrix, init_vec)}; 

		// 4. Assign new values to impurity.
		imp.set_vX(rotated_vec[0]);
		imp.set_vY(rotated_vec[1]);
		imp.set_vZ(rotated_vec[2]);
	}

	double calc_momentum_loss_freq(double ne, double te, 
		double n2, double t2, double v1, double m1, double m2, 
		double q1, double q2)
	{
		// This follows the derivation in Freidberg's textbook, Ch. 9. Most
		// everything is the same except the Coulomb logarithm is explicitly
		// calculated (i.e., the common assumption ln(alpha)=15 or 20 is not
		// enforced) and it has been generalized to any type of particles.
		// Generally use this as 1 = impurity, 2 = species colliding with. The
		// units of the input are:
		// ne, n2 [m-3]
		// te, t2 [eV]
		// v1 [m/s]
		// m1, m2 [kg]
		// q1, q2 [C]
		
		// Calculate reduced mass, Eq. 9.7
		double mr = calc_reduced_mass(m1, m2);

		// Debye length
		double lambda_debye {calc_debye_length(te, ne)};

		// The initial (squared) velocity of the reduced mass in the CoM 
		// frame. See Eq. 9.4 for the vector v. We square this for:
		// v^2 = v1^2 + v2^2 - 2 v1 dot v2
		// We assume that the second species is Maxwellian. So the last term
		// with the dot product vanishes on average since for every v2 vector
		// it is equally likely that an oppositely directed v2 vector exists.
		// So over many collisions, this term --> 0 on average. v is then
		// aligned along the x-axis in the CoM frame (Fig. 9.3) without any
		// loss in generality, thus v^2 --> v0^2.
		double mr_v0_sq {v1 * v1 + 1.5 * t2 
			* Constants::ev_to_j / m2};

		// b90, impact parameter that causes a 90 degree collision, Eq. 9.16
		double b90 {std::abs((q1 * q2) / (4.0 * Constants::pi 
			* Constants::eps0 * mr * mr_v0_sq))};

		// Coulomb factor, ln(alpha), see Eq. 9.34 but we retain some 
		// generality here instead of just approximating ln(alpha) = 20.
		double lnalpha {std::log(lambda_debye / b90)};

		// Thermal velocity (squared) of the species 1 is colliding 
		// with, Eq. 9.41. No 2/3 factor here because this comes after the 
		// problem has been simplified into a 1D problem, at least that's 
		// my understanding.
		double vt2_sq {2.0 * t2 * Constants::ev_to_j / m2};

		// The 1-2 collision frequency for the momentum loss rate. Compare to
		// Eq. 9.48. Only difference is e^4 --> q1^2 * q2^2.
		double nu {(1.0 / (4.0 * Constants::pi) * q1*q1 * q2*q2 * n2 / 
			(Constants::eps0*Constants::eps0 * m1 * mr) * lnalpha) *
			(1.0 / (v1*v1*v1 + 1.3 * vt2_sq))};

		return nu;
	}

	// Calculate a value for dt_coll, if needed. We use dt_coll as an
	// in/out variable here because it is not gauranteed that any 
	// restriction will be placed on the time step (e.g., if the particle is
	// at rest we don't have anything to say here).
	void set_var_time_step_coll(double& dt_coll, const Impurity::Impurity& imp,
		const Background::Background& bkg, const int tidx, const int xidx, 
		const int yidx, const int zidx, const Options::Options& opts)
	{
		// Calculate the momentum slowing down frequency. This would normally 
		// be done in collision_update, but we need to know the frequency to
		// calculate the time step so it gets a little mixed up.
	
		// Not handling neutral collisions (this is an approximation I do not
		// know if it is true or not, but it's common enough among codes). 
		if (imp.get_charge() == 0)
		{
			std::cerr << "Error! get_var_time_step_coll was called for a " 
				<< "neutral. Neutral-impurity collisions are not considered."
				<< " This is a programming error.\n";
				return;
		}

		// Calculate impurity velocity
		double imp_v {calc_imp_vel(imp)};

		// If impurity is at rest there's nothing to be done so return. We'd
		// get NaN's if we continued. Assign the minimum time step value.
		//if (imp_v < Constants::small) return;
		if (imp_v < Constants::small) dt_coll = opts.imp_time_step_min();

		// Load these into local variables since we pass them twice below.
		const double te = bkg.get_te()(tidx, xidx, yidx, zidx);
		const double ti = bkg.get_ti()(tidx, xidx, yidx, zidx);
		const double ne = bkg.get_ne()(tidx, xidx, yidx, zidx);

		// Load electron mass in kg. Calculate momentum loss frequency due to 
		// impurity-electron collisions. Note that we are passing in as the
		// correct units (J for temperatures, kg for masses and C for charges).
		double me {opts.gkyl_elec_mass_amu() * Constants::amu_to_kg};
		double nu_ze {calc_momentum_loss_freq(ne, te, ne, te, imp_v, 
			imp.get_mass(), me, -imp.get_charge() * Constants::charge_e, 
			Constants::charge_e)};

		// Load ion mass in kg. Calculate momentum loss frequency due to 
		// impurity-ion collisions. Assuming singly charged ions, but the 
		// function is generalized to any charge. Just need to change what's
		// passed in for q2 here if you want other charges. Also assuming
		// quasi-neutrality and passing in ni=ne here. 
		double mi {opts.gkyl_ion_mass_amu() * Constants::amu_to_kg};
		double nu_zi {calc_momentum_loss_freq(ne, te, ne, ti, imp_v, 
			imp.get_mass(),	mi, -imp.get_charge() * Constants::charge_e, 
			-Constants::charge_e)};

		// The change in velocity and the fraction are defined as: 
		//   imp_dv = -(nu_ze + nu_zi) * imp_v * imp_time_step
		//   mom_loss_frac = (imp_v + imp_dv) / imp_v
		// We don't want this changing by too much at one time, so we set a
		// reasonable limit on what the fraction can be and then calculate
		// the time step from that. Solve the above equations for dt to get:
		//   dt = (1 - f) / (nu_ze + nu_zi)
		double mom_loss_frac_limit {0.95};
		dt_coll = (1.0 - mom_loss_frac_limit) / (nu_ze + nu_zi);
		return;
	}

	void collision_update(Impurity::Impurity& imp, double te, double ti, 
		double ne, double imp_time_step, const Options::Options& opts,
		const bool split_particle, std::vector<Impurity::Impurity>& imps)
	{
		// Not handling neutral collisions (this is an approximation I do not
		// know if it is true or not, but it's common enough among codes). 
		if (imp.get_charge() == 0) return;

		// Calculate impurity velocity
		double imp_v {calc_imp_vel(imp)};

		// If impurity is at rest there's nothing to be done so return. We'd
		// get NaN's if we continued. 
		if (imp_v < Constants::small) return;

		// Load electron mass in kg. Calculate momentum loss frequency due to 
		// impurity-electron collisions. Note that we are passing in as the
		// correct units (J for temperatures, kg for masses and C for charges).
		double me {opts.gkyl_elec_mass_amu() * Constants::amu_to_kg};
		double nu_ze {calc_momentum_loss_freq(ne, te, ne, te, imp_v, 
			imp.get_mass(), me, -imp.get_charge() * Constants::charge_e, 
			Constants::charge_e)};
		//std::cout << std::scientific << "nu_ze = " << nu_ze << '\n';

		// Load ion mass in kg. Calculate momentum loss frequency due to 
		// impurity-ion collisions. Assuming singly charged ions, but the 
		// function is generalized to any charge. Just need to change what's
		// passed in for q2 here if you want other charges. Also assuming
		// quasi-neutrality and passing in ni=ne here. 
		double mi {opts.gkyl_ion_mass_amu() * Constants::amu_to_kg};
		double nu_zi {calc_momentum_loss_freq(ne, te, ne, ti, imp_v, 
			imp.get_mass(),	mi, -imp.get_charge() * Constants::charge_e, 
			-Constants::charge_e)};
		//std::cout << std::scientific << "nu_zi = " << nu_zi << '\n';

		// Test that the equation gives expected values for e-i (just to pick
		// and example collision type). 
		//double e_v {std::sqrt(3.0 * te * Constants::ev_to_j / me)};
		//double nu_ei_test {1.33e5 * ne * 1.0e-20 / 
		//	std::pow(ti * 1.0e-3, 1.5)};
		//double nu_ei {calc_momentum_loss_freq(ne, te, ne, ti, e_v, me, mi, 
		//	Constants::charge_e, -Constants::charge_e)};
		//std::cout << "e_v = " << e_v << '\n';
		//std::cout << std::scientific << "nu_ei = " << nu_ei << " " 
		//	<< nu_ei_test << "\n";

		// Calculate total momentum (i.e., change in velocity since mass is 
		// constant) loss from each type of collision.
		// d/dt(mz*vz) = -nu * mz * vz  --> dvz = -nu * vz * dt
		double imp_dv {-(nu_ze + nu_zi) * imp_v * imp_time_step};
		//std::cout << "imp_dv = " << imp_dv << '\n';

		double mom_loss_frac {(imp_v + imp_dv) / imp_v};
		//std::cout << "mom_loss_frac = " << mom_loss_frac << '\n';
	
		//std::cout << "Before: " << imp.get_vx() << ", " << imp.get_vy() 
		//	<< ", " << imp.get_vz() << '\n';
		if (opts.imp_collisions_int() == 2)
		{
			// Elastic approach. Straightforward approach to apply this to the 
			// impurity is to scale down each impurity velocity component 
			// proportional to the amount that was lost in the collisions.
			// Before: v0 = (vx0, vy0, vz0)
			// After:  v1 = a * v0 = (a*vx0, a*vy0, a*vz0) = (vx1, vy1, vz1)
			// Thus each component of the after-collision vector is:
			// vx1 = a * vx0    and likewise for y and z.
			// In theory this can be negative if dv < v. This would just mean the
			// particle is turning around due to a collision I suppose. 
			imp.set_vX(imp.get_vX() * mom_loss_frac); 
			imp.set_vY(imp.get_vY() * mom_loss_frac); 
			imp.set_vZ(imp.get_vZ() * mom_loss_frac); 
		}
		else if (opts.imp_collisions_int() == 1)
		{
			// Inelastic approach. Consider a new coordinate system rotated by
			// by some unknown angle such that the particle velocity is along
			// the x axis (this assumption was actually already made in
			// in calculating the momentum loss frequency). 
			// Then before/after a collision looks like:
			// Before: v0 = (v0x,   0,   0)
			// After:  v1 = (v1x, v1y, v1z)
			// The calculated imp_dv is the velocity (really momentum, but 
			// mass is constant so we can just talk velocities) loss over this 
			// time step **in the x direction**. So we can rewrite v1 as:
			//   v1 = (v0x + dv, v1y, v1z)
			// The angle between the before/after vectors is the dot product:
			//   v0 * v1 = |v0||v1|cos(theta) = v0x * v1x  + 0 + 0
			// We then assume an inelastic collision approximation,  
			// so |v0||v1| = v0x^2. This is because we do not have enough
			// information to know how much |v1| changes by, we only know how
			// much the x component changed. Thus:
			//   cos(theta) = v0x * v1x / v0x^2 = v1x / v0x 
			//     = (v0x + dv) / v0x = mom_loss_frac
			// So we can achieve the post-collision vector by just rotating the
			// original vector by theta. The next question is then in which
			// direction - x, y, z or some combination of the three? Well this
			// is a random process in which *any* direction is equally
			// possible, so we simply create a random unit vector and rotate
			// our vector around it by theta. This is to say, we assume the 
			// species it is colliding with is Maxwellian, which is fine
			// enough. Once we have theta (defl_ang), we rotate the particle
			// vector by it and we're done. 
			// There is a subtlety in that theta was calculated after we had
			// chosen a coordinate system aligned with the x direction. This
			// is fine, because we just need to know what angle to deflect by,
			// which will be the same angle in both the original or our 
			// temporary coordinate system (theta can be thought of as an 
			// displacement - if you're moving 5 meters, it doesn't matter if
			// you did it in Virginia or California, it's the same either way).
			double defl_ang {std::acos(mom_loss_frac)};

			// Equal chance of being positive or negative angle. No point in
			// doing this since the particle is rotated around a random axis.
			//if (Random::get(0, 1) == 1) defl_ang *= -1.0;

			// If particle is to be split, then evenly split the particle's 
			// weight with a secondary particle and have that particle rotated 
			// with the opposite deflection angle. Evenly split because there's
			// equal chance of either deflection angle.
			if (split_particle)
			{
				// Split particle is appended to end of imps to be followed
				// later. 
				VarianceReduction::create_secondary(imp, imps, 
					imp.get_weight() / 2.0);

				// Don't forget to rotate the secondary (back() returns
				// a reference to the split particle that was just appended
				// onto the end of imps in create_secondary)
				rotate_imp(imps.back(), defl_ang);
			}

			//std::cout << "defl_ang = " << defl_ang << '\n';
			rotate_imp(imp, defl_ang);
		}
		//std::cout << "After:  " << imp.get_vx() << ", " << imp.get_vy() 
		//	<< ", " << imp.get_vz() << '\n';
	}

	std::tuple<double, double, double> sample_bkg_velocity(
		const Background::Background& bkg, const int tidx, const int xidx, 
		const int yidx, const int zidx, const double T, const double m)
	{
		// Get mean parallel flow, which is the drifting term in a drifting
		// Maxwellian. This is the ion velocity, which assume is the plasma
		// velocity (so same for both ions and electrons). 
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

	std::tuple<double, double, double, double> nanbu_calc_s(
		Impurity::Impurity& imp, const Background::Background& bkg, 
		const int tidx, const int xidx, const int yidx, const int zidx, 
		const double imp_time_step, const double T, const double mass_kg, 
		const double Te, const double ne)
	{
		// Take random sample of background species being collided with
		auto [bkg_vX, bkg_vY, bkg_vZ] = sample_bkg_velocity(bkg, tidx, xidx, 
			yidx, zidx, T, mass_kg);
		//std::cout << "T = " << T << '\n';
		//std::cout << "Te = " << Te << '\n';
		//std::cout << "ne = " << ne << '\n';
		//std::cout << "bkg_vX = " << bkg_vX << '\n';
		//std::cout << "bkg_vY = " << bkg_vY << '\n';
		//std::cout << "bkg_vZ = " << bkg_vZ << '\n';
		//std::cout << "imp_vX = " << imp.get_vX() << '\n';
		//std::cout << "imp_vY = " << imp.get_vY() << '\n';
		//std::cout << "imp_vZ = " << imp.get_vZ() << '\n';

		// Relative velocity, g
		double gX {imp.get_vX() - bkg_vX};
		double gY {imp.get_vY() - bkg_vY};
		double gZ {imp.get_vZ() - bkg_vZ};
		double g {std::sqrt(gX*gX + gY*gY + gZ*gZ)};
		//std::cout << "g = " << g << '\n';
		//if (g < Constants::small) g = Constants::small;
		if (g < Constants::small) g = 0.001;

		// Reduced mass of impurity and species colliding with
		double mu_ab {imp.get_mass() * mass_kg / (imp.get_mass() + mass_kg)};
		//std::cout << "mu_ab = " << mu_ab << '\n';

		// Expectation value of g^2
		double expect_g_sq {3.0 * T * Constants::ev_to_j / mu_ab};
		//std::cout << "expect_g_sq = " << expect_g_sq << '\n';
		//if (expect_g_sq < Constants::small) expect_g_sq = Constants::small;
		if (expect_g_sq < Constants::small) expect_g_sq = 0.001;

		// Expectation value of b0. Assuming singly charged background species.
		double expect_b0 {imp.get_charge() * Constants::charge_e 
			* Constants::charge_e / (2.0 * Constants::pi * Constants::eps0 
			* mu_ab * expect_g_sq)};
		//std::cout << "expect_b0 = " << expect_b0 << '\n';

		// Debye length and Coulumb logarithm
		double debye_length {calc_debye_length(Te, ne)};
		double ln_alpha {std::log(debye_length / expect_b0)};
		if (std::isnan(ln_alpha))
		{
			std::cerr << "Error! ln_alpha = nan\n";
			std::cerr << "Te = " << Te << '\n';
			std::cerr << "ne = " << ne << '\n';
			std::cerr << "debye_length = " << debye_length << '\n';
			std::cerr << "expect_b0 = " << expect_b0 << '\n';
			std::cerr << "ln_alpha = " << ln_alpha << '\n';
		}

		// Calculate s (Eq. 19). Assume ni=ne and singly charge background
		double square_term {imp.get_charge() * Constants::charge_e 
			* Constants::charge_e / (Constants::eps0 * mu_ab)};
		//std::cout << "square_term = " << square_term << '\n';
		double s {ln_alpha / (4.0 * Constants::pi) * square_term * square_term
			* ne / (g*g*g) * imp_time_step};

		if (std::isnan(s))
		{
			std::cerr << "Error! s = nan\n";
			std::cerr << "T = " << T << '\n';
			std::cerr << "Te = " << Te << '\n';
			std::cerr << "bkg_vX = " << bkg_vX << '\n';
			std::cerr << "bkg_vY = " << bkg_vY << '\n';
			std::cerr << "bkg_vZ = " << bkg_vZ << '\n';
			std::cerr << "imp_vX = " << imp.get_vX() << '\n';
			std::cerr << "imp_vY = " << imp.get_vY() << '\n';
			std::cerr << "imp_vZ = " << imp.get_vZ() << '\n';
			std::cerr << "g = " << g << '\n';
			std::cerr << "expect_g_sq = " << expect_g_sq << '\n';
			std::cerr << "ne = " << ne << '\n';
			std::cerr << "debye_length = " << debye_length << '\n';
			std::cerr << "expect_b0 = " << expect_b0 << '\n';
			std::cerr << "ln_alpha = " << ln_alpha << '\n';

		}

		return std::make_tuple(s, gX, gY, gZ);
	}

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
			// calculated using python/calc_nabu.py.
			return Utilities::linear_interpolate(Nanbu::s, Nanbu::A, s);
		}

	}

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
		//std::cout << "cos_chi = " << cos_chi << '\n';
		return std::acos(cos_chi);

	}

	std::tuple<double, double, double> nanbu_post_coll(const double gX,
		const double gY, const double gZ, const double chi, const double mass_a,
		const double mass_b, const Impurity::Impurity& imp)
	{
		// For calculating h components
		double g_perp {std::sqrt(gY*gY + gZ*gZ)};
		double g {std::sqrt(gX*gX + gY*gY + gZ*gZ)};
		double eps {2.0 * Constants::pi * Random::get(0.0, 1.0)};

		// Calculate the h components
		double hX {g_perp * std::cos(eps)};
		double hY {-(gY * gX * std::cos(eps) 
			+ g * gZ * std::sin(eps)) / g_perp};
		double hZ {-(gZ * gX * std::cos(eps) 
			- g * gY * std::sin(eps)) / g_perp};

		// Calculate post collision velocities
		double mu {mass_b / (mass_a + mass_b)};
		double cos_chi {std::cos(chi)};
		double sin_chi {std::sin(chi)};
		double vX_post {imp.get_vX() - mu * (gX 
			* (1.0 - cos_chi) + hX * sin_chi)};
		double vY_post {imp.get_vY() - mu * (gY 
			* (1.0 - cos_chi) + hY * sin_chi)};
		double vZ_post {imp.get_vZ() - mu * (gZ 
			* (1.0 - cos_chi) + hZ * sin_chi)};

		return std::make_tuple(vX_post, vY_post, vZ_post);
	}

	void nanbu_coll(Impurity::Impurity& imp, const Background::Background& bkg,
		const int tidx, const int xidx, const int yidx, const int zidx,
		const Options::Options& opts, bool elec, const double imp_time_step)
	{
		// Derivation and steps taken from:
		// Nanbu, K. Theory of cumulative small-angle collisions in plasmas. 
		// Phys. Rev. E 55, 4642–4652 (1997).

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
		//std::cout << "s = " << s << '\n';
		//std::cout << "gX = " << gX << '\n';
		//std::cout << "gY = " << gY << '\n';
		//std::cout << "gZ = " << gZ << '\n';

		// Calcluate A (Eq. 13)
		double A {nanbu_calc_A(s)};
		//std::cout << "A = " << A << '\n';

		// Calculate deflection angle, chi (Eq. 17)
		double chi {nanbu_calc_chi(s, A)};
		if (std::isnan(chi))
		{
			std::cerr << "Error! chi = nan\n";
			std::cerr << "T = " << T << '\n';
			std::cerr << "Te = " << Te << '\n';
			std::cerr << "ne = " << ne << '\n';
			std::cerr << "chi = " << chi << '\n';
			std::cerr << "s = " << s << '\n';
			std::cerr << "A = " << A << '\n';
		}

		// Calculate post-collision velocity (Eq. 20a)
		auto [vX_post, vY_post, vZ_post] = nanbu_post_coll(gX, gY, gZ, chi, 
			imp.get_mass(), mass_kg, imp);
		//std::cout << "pre: vX, vY, vZ = " << imp.get_vX() << ", " 
		//	<< imp.get_vY() << ", " << imp.get_vZ() << '\n';
		//std::cout << "post: vX, vY, vZ = " << vX_post << ", " << vY_post << ", " 
		//	<< vZ_post << '\n';

		// Error if velocity magnitude changes by more than X%
		/*
		double pre_imp_v {std::sqrt(imp.get_vX() * imp.get_vX() + 
			imp.get_vY() * imp.get_vY() + imp.get_vZ() * imp.get_vZ())};
		double post_imp_v {std::sqrt(vX_post * vX_post + vY_post * vY_post 
			+ vZ_post * vZ_post)};
		if (std::abs((post_imp_v - pre_imp_v) / pre_imp_v) > 0.1)
		{
			std::cerr << "Large collision detected!\n";
			std::cerr << "pre: vX, vY, vZ = " << imp.get_vX() << ", " 
				<< imp.get_vY() << ", " << imp.get_vZ() << '\n';
			std::cerr << "post: vX, vY, vZ = " << vX_post << ", " << vY_post << ", " 
				<< vZ_post << '\n';
			std::cerr << "T = " << T << '\n';
			std::cerr << "Te = " << Te << '\n';
			std::cerr << "ne = " << ne << '\n';
			std::cerr << "s = " << s << '\n';
			std::cerr << "A = " << A << '\n';
			std::cerr << "chi = " << chi << '\n';
		}
		*/

		// Update impurity
		imp.set_vX(vX_post);
		imp.set_vY(vY_post);
		imp.set_vZ(vZ_post);


	}
}
