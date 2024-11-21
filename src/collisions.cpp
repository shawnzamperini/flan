/**
* @file collisions.cpp
*/
#include <cmath>
#include <utility>
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>

#include "impurity.h"
#include "read_input.h"
#include "constants.h"
#include "random.h"

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

	double calc_imp_vel(Impurity::Impurity& imp)
	{
		return std::sqrt(imp.get_vx() * imp.get_vx() 
			+ imp.get_vy() * imp.get_vy() + imp.get_vz() * imp.get_vz());
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
		// to randians first.
		std::array<std::array<double, 3>, 3> rotation_matrix {
			calc_rotation_matrix(rotation_axis, 
			defl_ang)};

		// 3. Apply rotation matrix to find new rotated vector
		std::array<double, 3> init_vec {imp.get_vx(), imp.get_vy(), 
			imp.get_vz()};
		std::array<double, 3> rotated_vec {multiply_matrix_vector(
			rotation_matrix, init_vec)}; 

		// 4. Assign new values to impurity.
		imp.set_vx(rotated_vec[0]);
		imp.set_vy(rotated_vec[1]);
		imp.set_vz(rotated_vec[2]);
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

	void collision_update(Impurity::Impurity& imp, double te, double ti, 
		double ne, double imp_time_step)
	{
		// Not handling neutral collisions (this is an approximation I do not
		// know if it is true or not, but it's common enough among codes). 
		if (imp.get_charge() == 0) return;

		//std::cout << "=== Collisions ===\n";

		// Calculate impurity velocity
		double imp_v {calc_imp_vel(imp)};
		//std::cout << "imp_v = " << imp_v << '\n';

		// If impurity is at rest there's nothing to be done so return. We'd
		// get NaN's if we continued. 
		if (imp_v < Constants::small) return;

		// Load electron mass in kg. Calculate momentum loss frequency due to 
		// impurity-electron collisions. Note that we are passing in as the
		// correct units (J for temperatures, kg for masses and C for charges).
		double me {Input::get_opt_dbl(Input::gkyl_elec_mass_amu) 
			* Constants::amu_to_kg};
		double nu_ze {calc_momentum_loss_freq(ne, te, ne, te, imp_v, 
			imp.get_mass(), me, -imp.get_charge() * Constants::charge_e, 
			Constants::charge_e)};
		//std::cout << std::scientific << "nu_ze = " << nu_ze << '\n';

		// Load ion mass in kg. Calculate momentum loss frequency due to 
		// impurity-ion collisions. Assuming singly charged ions, but the 
		// function is generalized to any charge. Just need to change what's
		// passed in for q2 here if you want other charges. Also assuming
		// quasi-neutrality and passing in ni=ne here. 
		double mi {Input::get_opt_dbl(Input::gkyl_ion_mass_amu) 
			* Constants::amu_to_kg};
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

		// Straightforward approach to apply this to the impurity is to scale 
		// down each impurity velocity component proportional to the amount
		// that was lost in the collisions.
		// Before: v0 = (vx0, vy0, vz0)
		// After:  v1 = a * v0 = (a*vx0, a*vy0, a*vz0) = (vx1, vy1, vz1)
		// Thus each component of the after-collision vector is:
		// vx1 = a * vx0    and likewise for y and z.
		// In theory this can be negative if dv < v. This would just mean the
		// particle is turning around due to a collision I suppose. 
		double mom_loss_frac {(imp_v + imp_dv) / imp_v};
	
		// Just a temporary way to choose between inelastic (true) or elastic
		// (false) collisions.
		if (true)
		{
			std::cout << "mom_loss_frac = " << mom_loss_frac << '\n';
			//std::cout << "Before: " << imp.get_vx() << ", " << imp.get_vy() 
			//	<< ", " << imp.get_vz() << '\n';
			imp.set_vx(imp.get_vx() * mom_loss_frac); 
			imp.set_vy(imp.get_vy() * mom_loss_frac); 
			imp.set_vz(imp.get_vz() * mom_loss_frac); 
			//std::cout << "After:  " << imp.get_vx() << ", " << imp.get_vy() 
			//	<< ", " << imp.get_vz() << '\n';
		}

		// This gives weird results, probably remove eventually. 
		else
		{
			// Inelastic approach. The axes are rotated by some unknown angle 
			// so 
			// that the initial velocity lies along x. This was already done/
			// impied in the previous steps, just stating here that it happened.
			// Then before/after looks like:
			// Before: v0 = (v0x,   0,   0)
			// After:  v1 = (v1x, v1y, v1z)
			// The calculated imp_dv is the velocity (really momentum, but 
			// mass is
			// constant so we can just talk velocities) loss over this time 
			// step
			// **in the direction of the vector**. So we can rewrite v1 as:
			//   v1 = (v0x - dv, v1y, v1z)
			// The angle between the vector is the dot product:
			//   v0 * v1 = |v0||v1|cos(theta) = v0x * v1x 
			// This is an inelastic collision approximation, so |v0| = |v1| 
			// = v0x^2
			// thus:
			//   cos(theta) = v1x / v0x = (v0x - dv) / v0x
			// So we can achieve the post-collision vector by just rotating the
			// original vector by theta. But what in what direction? And what 
			// about
			// the original rotation that we did? Ignore the original because 
			// fuck
			// it, end of the day we just need to rotate the vector in a random
			// direction by theta since any direction is equally possible, 
			// assuming
			// the species the impurity is colliding with is Maxwellian (which we
			// do assume here).
			//std::cout << "Before: " << imp.get_vx() << ", " << imp.get_vy() 
			//	<< ", " << imp.get_vz() << '\n';
			double defl_ang {std::acos(mom_loss_frac)};
			//std::cout << "defl_ang = " << defl_ang << '\n';
			rotate_imp(imp, defl_ang);
			//std::cout << "After:  " << imp.get_vx() << ", " << imp.get_vy() 
			//	<< ", " << imp.get_vz() << '\n';
		}
	}
}
