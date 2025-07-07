/**
* @file collisions.cpp
*/
#include <cmath>
#include <utility>
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>

#include "background.h"
#include "constants.h"
#include "impurity.h"
#include "options.h"
#include "random.h"
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
}
