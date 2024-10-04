/**
* @file collisions.cpp
*/
#include <cmath>
#include <utility>
#include <vector>
#include <array>
#include <iostream>

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
		return std::sqrt(Constants::eps0 * te / 
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
		}

		return vels;
	}

	double calc_imp_vel(Impurity::Impurity& imp)
	{
		return std::sqrt(imp.get_vx() * imp.get_vx() 
			+ imp.get_vy() * imp.get_vy() + imp.get_vz() * imp.get_vz());

	}
	
	std::vector<double> generate_impact_param(double bmin, double bmax,
		int num)
	{
		std::vector<double> b (num);
		for (int i {}; i < num; ++i)
		{
			// This is the cumulative distribution function for the pdf 1/r.
			b[i] = bmax - bmin * std::exp(Random::get(0.0, 1.0) 
				* std::log(bmax / bmin));
		}
		return b;
	}

	double calc_defl_ang(double red_mass, std::vector<double>& vel, 
		std::vector<double>& b, int q1, int q2)
	{
		// Calculate deflection angle over the number of collisions specified
		// by the length of vel (or b, same size). 
		double defl_ang {};
		for (std::size_t i {}; i < vel.size(); ++i)
		{
			defl_ang += 2.0 * std::atan(1.0 / (4.0 * Constants::pi 
				* Constants::eps0 * red_mass / (q1 * q2) * vel[i] * vel[i] 
				* b[i]));
		}
		return defl_ang;
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
			defl_ang * Constants::pi / 180.0)};

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

	void collision_step(Impurity::Impurity& imp, 
		const double te, const double ti, const double ne, 
		double imp_time_step)
	{

		std::cout << "=============================\n";
		// Load some input options up front
		double me {Input::get_opt_dbl(Input::gkyl_elec_mass_amu) 
			* Constants::amu_to_kg};
		double mi {Input::get_opt_dbl(Input::gkyl_ion_mass_amu)
			* Constants::amu_to_kg};

		// First get the reduced masses for each collision option. Right now 
		// this is specific to collisions with electrons and an ion species
		double red_mass_ze = calc_reduced_mass(imp.get_mass(), me);
		double red_mass_zi = calc_reduced_mass(imp.get_mass(), mi);

		// Then calculate the Debye length. Only collisions within a radius 
		// of this are considered.
		double lambda_debye {calc_debye_length(te, ne)};
		std::cout << "lambda_debye = " << lambda_debye << '\n';

		// Sample a large number of different electron and ion velocites from
		// a Maxwellian distribution centered on Te, Ti. The algoritm returns
		// two numbers at a time.
		int imp_coll_num {Input::get_opt_int(Input::imp_coll_num)};

		// Calculate average ion velocity. We don't really have the standard 
		// deviation of the the ion velocity, so we're hardcoding it as avg/4.
		// This could certainly be improved in the future. Multiplying by
		// ev_to_j converts eV --> Joules.
		double avg_ion_vel {2.0 * ti * Constants::ev_to_j / mi};
		double std_ion_vel {avg_ion_vel / 4.0};
		std::vector<double> ion_vels {generate_random_vels(avg_ion_vel, 
			std_ion_vel, imp_coll_num)};
		std::cout << "avg_ion_vel = " << avg_ion_vel << '\n';

		// Likewise for electron velocity.
		double avg_elec_vel {2.0 * te * Constants::ev_to_j / me};
		double std_elec_vel {avg_elec_vel / 4.0};
		std::vector<double> elec_vels {generate_random_vels(avg_elec_vel, 
			std_elec_vel, imp_coll_num)};
		std::cout << "avg_elec_vel = " << avg_elec_vel << '\n';

		// Then calculate the corresponding center of mass velocities. This
		// needs the impurity velocity. We are reusing ion_vels and elec_vels
		// to save some time.
		double vz {calc_imp_vel(imp)};
		for (std::size_t i {}; i < ion_vels.size(); ++i)
		{
			ion_vels[i] = std::sqrt(vz + ion_vels[i] * ion_vels[i]);
		}
		for (std::size_t i {}; i < elec_vels.size(); ++i)
		{
			elec_vels[i] = std::sqrt(vz + elec_vels[i] * elec_vels[i]);
		}

		// Genarate a random number of impact parameters (b) from a 1/r
		// distribution between some specified bmin and the Debye length.
		// bmin just needs to be artbitrarily small, so right now it is 
		// hardcoded as 1e-10 m. 
		double bmin {1e-10};
		std::vector<double> ion_b {generate_impact_param(bmin, 
			lambda_debye, imp_coll_num)}; 
		std::vector<double> elec_b {generate_impact_param(bmin, 
			lambda_debye, imp_coll_num)}; 
		
		// Now calculate a vector containing all the deflection angles. The
		// highest fidelity approach would be to see how many total ion/
		// electron pairs the impurity would encounter and then calculate that
		// many mini-deflection angles. This would just take too long, so 
		// instead we calculate a smaller number of deflection (imp_coll_num),
		// enough to generate good statistics, and then artifically scale it
		// up to the total number of particles by multiplying to the total
		// deflection angle after imp_coll_num collisions by real_num_coll /
		// imp_coll_num. 

		// Deflection angle for elec-imp collisions
		double elec_imp_ang {calc_defl_ang(red_mass_ze, elec_vels, elec_b, 
			Constants::charge_e, -imp.get_charge() * Constants::charge_e)};
		std::cout << "elec_imp_ang = " << elec_imp_ang << '\n';

		// Deflection angle for ion-imp collisions. Assuming singly charged
		// ions here, though the function is general and can handle 
		// other charge states if that's needed.
		double ion_imp_ang {calc_defl_ang(red_mass_zi, ion_vels, ion_b, 
			-Constants::charge_e, -imp.get_charge() * Constants::charge_e)};
		std::cout << "ion_imp_ang = " << ion_imp_ang << '\n';

		// At this point we have a deflection angle for an arbitrary 
		// imp_coll_num number of collisions. Over this timestep, the impurity
		// likely encounters many more collisions though. So to estimate it
		// without setting imp_coll_num to crazy large numbers, calculate by
		// what factor the real number of collisions is larger by, and multiply
		// the angle we have by that. 
	
		// Calculate how many electron/ion pairs the impurity "sees" in the
		// Debye cylinder corresponding to the impurity's distance traveled
		// during the time step (dist = vz * imp_time_step). 
		double n_debye_cyl {Constants::pi * Constants::eps0 * te 
			* imp_time_step * vz / (Constants::charge_e * Constants::charge_e)};
		std::cout << "n_debye_cyl = " << n_debye_cyl << '\n';

		double tot_defl_ang {(n_debye_cyl / static_cast<double>(imp_coll_num)) 
			* (elec_imp_ang + ion_imp_ang)};	

		// Now modify the impurity ion's trajectory by tot_defl_ang.
		std::cout << "imp velocity\n";
		std::cout << imp.get_vx() << ", " << imp.get_vy() << ", " 
			<< imp.get_vz() << '\n';
		rotate_imp(imp, tot_defl_ang);
		std::cout << imp.get_vx() << ", " << imp.get_vy() << ", " 
			<< imp.get_vz() << '\n';

	}
}
