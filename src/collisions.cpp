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
	
	std::vector<double> generate_impact_param(double bmin, double bmax,
		int num)
	{
		std::vector<double> b (num);
		for (int i {}; i < num; ++i)
		{
			// This is the cumulative distribution function for the pdf 1/r.
			//b[i] = bmax - bmin * std::exp(Random::get(0.0, 1.0) 
			//	* std::log(bmax / bmin));

			// Uniform distribution function.
			//b[i] = bmin + Random::get(0.0, 1.0) * (bmax - bmin);

			// Cumulative distribution function for a 1/r^2 pdf
			b[i] = std::sqrt(1.0 / (Random::get(0.0, 1.0) 
				* (1.0 / bmin - 1.0 / bmax)));
		}
		return b;
	}

	double calc_defl_freq(double red_mass, std::vector<double>& vel, 
		std::vector<double>& b, double q1, double q2, double lambda_debye,
		double vz, double imp_time_step, double te, double ne)
	{
		// Calculate deflection angle over the number of collisions specified
		// by the length of vel (or b, same size). 
		double defl_ang {};
		//std::cout << std::scientific << "  q1 = " << q1 << '\n';
		//std::cout << std::scientific << "  q2 = " << q2 << '\n';
		//std::cout << "  red_mass = " << red_mass << '\n';
		for (std::size_t i {}; i < vel.size(); ++i)
		{
			//std::cout << "  b[" << i << "] = " << b[i] << '\n';
			//std::cout << "  vel[" << i << "] = " << vel[i] << '\n';
			//std::cout << "  defl_ang term = " << 2.0 * 1.0 / 
			//	(4.0 * Constants::pi 
			//	* Constants::eps0 * red_mass / (q1 * q2) * vel[i] * vel[i] 
			//	* b[i]) << '\n';
			//std::cout << "  defl_ang = " << 2.0 * std::atan(1.0 / 
			//	(4.0 * Constants::pi 
			//	* Constants::eps0 * red_mass / (q1 * q2) * vel[i] * vel[i] 
			//	* b[i])) << '\n';
			defl_ang += 2.0 * std::atan(1.0 / (4.0 * Constants::pi 
				* Constants::eps0 * red_mass / (q1 * q2) * vel[i] * vel[i] 
				* b[i]));
		}

		// The timescale over which these collisions occur is how long it takes
		// the impurity to leave the Debye sphere, i.e., dt = (lambda / 2) / vz.
		double dt {lambda_debye / 2.0 / vz};

		// Calculate average deflection frequency
		double avg_defl_per_coll {defl_ang / static_cast<double>(vel.size())};
		std::cout << "method #1: " << avg_defl_per_coll / dt << '\n';
		//double avg_defl_freq {avg_defl_per_coll / dt};
		//double n_debye_cyl {Constants::pi * Constants::eps0 * te 
		//	* Constants::ev_to_j * imp_time_step * vz / (Constants::charge_e 
		//	* Constants::charge_e)};
		double n_debye_cyl {Constants::pi * lambda_debye * lambda_debye * vz
			* imp_time_step * ne};
		double avg_defl_freq {defl_ang / static_cast<double>(vel.size())
			* n_debye_cyl / imp_time_step};
		std::cout << "method #2: " << avg_defl_per_coll / imp_time_step 
			* n_debye_cyl << '\n';

		return avg_defl_freq;
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

	std::vector<double> generate_com_vels(Impurity::Impurity& imp, 
		std::vector<double>& vels)
	{
		// Vector to hold CoM velocities.
		std::vector<double> com_vels (vels.size());

		for (std::size_t i {}; i < vels.size(); ++i)
		{
			// Each loop represents an ion/electron that we need to generate
			// a randomly oriented vector with magnitude vel for. So get a 
			// random unit vector first.
			std::array<double, 3> unit_vector {ran_unit_vector()};

			// The magnitude of the CoM velocity vector is then
			com_vels[i] = std::sqrt(
			std::pow(imp.get_vx() - unit_vector[0] * vels[i], 2)
			+ std::pow(imp.get_vy() - unit_vector[1] * vels[i], 2)
			+ std::pow(imp.get_vz() - unit_vector[2] * vels[i], 2));
		}
		return com_vels;
	}

	void compare_theory_90_zi(double calc_freq, double ne, double te, 
		double ti, double vz, double red_mass, double mz, double mi, 
		double qz, double qi)
	{
		// This helper function compares the given angular deflection frequency
		// to the 90-degree frequency predicted by theory for ion-impurity
		// collisions. It is just to tell us if what we are doing makes sense 
		// in regards to established theory.
		// ne [m-3]
		// te, ti [K]
		// vz [m/s]
		// red_mass [kg]
		// mz, mi [kg]
		// qz, qi [C]

		// Calculate ln(alpha)
		double lnalpha_fact {std::log(std::pow(Constants::eps0, 1.5) * 4.0 *
			Constants::pi * Constants::amu_to_kg / 
			std::pow(Constants::ev_to_j, 2.5))};
		double lnalpha {lnalpha_fact + std::log(std::sqrt(te / 
			Constants::ev_to_j / ne) * vz * vz * (red_mass 
			/ Constants::amu_to_kg) / (std::abs(qz / Constants::ev_to_j 
			* qi / Constants::ev_to_j)))};

		// Calculate the 90 degree deflection frequency
		double nu90_fact {std::pow(Constants::ev_to_j, 4.0) * 
			std::pow(Constants::amu_to_kg, 1.5) / (4.0 * Constants::pi 
			* std::pow(Constants::eps0, 2.0) 
			* std::pow(Constants::amu_to_kg, 2.0) 
			* std::pow(2.0 * Constants::ev_to_j, 1.5))};
		double nu90 {nu90_fact * std::pow(qz / Constants::ev_to_j, 2.0) 
		* std::pow(qi / Constants::ev_to_j, 2.0) * ne * lnalpha / 
		((mz / Constants::amu_to_kg) * (red_mass / Constants::amu_to_kg))
		/ std::pow((ti / Constants::ev_to_j) / (mz / Constants::amu_to_kg), 
		1.5) + 1.3 * std::pow((ti / Constants::ev_to_j) / 
		(mi / Constants::amu_to_kg), 1.5)};

		// Calculate what the 90 degree deflection frequency is for our case
		double calc_nu90 {calc_freq / (Constants::pi / 2.0)};

		// Print out
		std::cout << std::scientific << "theory, calc = " << nu90 
			<< ", " << calc_nu90 << '\n';

	}

	double calc_momentum_loss_freq(double ne, double te, 
		double n2, double t2, double v1, double m1, double m2, 
		double q1, double q2)
	{
		// Generally 1 = impurity, 2 = species colliding with
		// ne [m-3]
		// te, temp2 [eV]
		// v1 [m/s]
		// m1, m2 [kg]
		// q1, q2 [C]
		
		// Calculate reduced mass
		double mr = calc_reduced_mass(m1, m2);

		// Debye length
		double lambda_debye {calc_debye_length(te, ne)};

		// The initial (squared) velocity of the reduced mass in the CoM 
		// frame. See derivation in ...
		double mr_v0_sq {v1 * v1 + 1.5 * t2 
			* Constants::ev_to_j / m2};

		// b90, impact parameter that causes a 90 degree collision
		double b90 {std::abs((q1 * q2) / (4.0 * Constants::pi 
			* Constants::eps0 * mr * mr_v0_sq))};

		// Coulomb factor, ln(alpha)
		double lnalpha {std::log(lambda_debye / b90)};

		// Thermal velocity (squared) of the species 1 is colliding with. No
		// 2/3 factor here because this comes after the problem has been 
		// simplified into a 1D problem, at least that's my understanding.
		double vt2_sq {2.0 * t2 * Constants::ev_to_j / m2};

		// The 1-2 collision frequency for the momentum loss rate.
		double nu {(1.0 / (4.0 * Constants::pi) * q1*q1 * q2*q2 * n2 / 
			(Constants::eps0*Constants::eps0 * m1 * mr) * lnalpha) *
			(1.0 / (v1*v1*v1 + 1.3 * vt2_sq))};

		return nu;
	}

	double old_calc_momentum_loss_freq(double ne, double te, 
		double temp2, double v1, double m1, double m2, 
		double q1, double q2)
	{
		// Generally 1 = impurity, 2 = species colliding with
		// ne [m-3]
		// te, temp2 [J]
		// v1 [m/s]
		// red_mass [kg]
		// m1, m2 [kg]
		// q1, q2 [C]
		
		// Calculate reduced mass
		double red_mass = calc_reduced_mass(m1, m2);

		// Calculate ln(alpha). Broken up to make it a bit more digestible. 
		double lnalpha_fact {std::log(std::pow(Constants::eps0, 1.5) * 4.0 *
			Constants::pi * Constants::amu_to_kg / 
			std::pow(Constants::ev_to_j, 2.5))};
		double lnalpha {lnalpha_fact + std::log(std::sqrt(te / 
			Constants::ev_to_j / ne) * v1 * v1 * (red_mass 
			/ Constants::amu_to_kg) / (std::abs(q1 / Constants::ev_to_j 
			* q2 / Constants::ev_to_j)))};
		std::cout << "  lnalpha_fact = " << lnalpha_fact << '\n';
		std::cout << "  lnalpha = " << lnalpha << '\n';

		// Calculate the momentum loss frequency. This is the frequency at
		// which the species 1 loses its forward-directed momentum due to
		// collisions with species 2.
		double nu_fact {std::pow(Constants::ev_to_j, 4.0) * 
			std::pow(Constants::amu_to_kg, 1.5) / (4.0 * Constants::pi 
			* std::pow(Constants::eps0, 2.0) 
			* std::pow(Constants::amu_to_kg, 2.0) 
			* std::pow(2.0 * Constants::ev_to_j, 1.5))};
		double nu {nu_fact * std::pow(q1 / Constants::ev_to_j, 2.0) 
		* std::pow(q2 / Constants::ev_to_j, 2.0) * ne * lnalpha / 
		((m1 / Constants::amu_to_kg) * (red_mass / Constants::amu_to_kg))
		/ std::pow((temp2 / Constants::ev_to_j) / (m1 / Constants::amu_to_kg), 
		1.5) + 1.3 * std::pow((temp2 / Constants::ev_to_j) / 
		(m2 / Constants::amu_to_kg), 1.5)};
		std::cout << "  nu_fact = " << nu_fact << '\n';
		std::cout << "  nu = " << nu << '\n';

		return nu;

	}

	void collision_step(Impurity::Impurity& imp, double te, double ti, 
		double ne, double imp_time_step)
	{
		
		std::cout << "=== Collisions ===\n";

		// Calculate impurity velocity
		double imp_v {calc_imp_vel(imp)};
		std::cout << "imp_v = " << imp_v << '\n';

		// Load electron mass in kg. Calculate momentum loss frequency due to 
		// impurity-electron collisions. Note that we are passing in as the
		// correct units (J for temperatures, kg for masses and C for charges).
		double me {Input::get_opt_dbl(Input::gkyl_elec_mass_amu) 
			* Constants::amu_to_kg};
		//double nu_ze {calc_momentum_loss_freq(ne, te * Constants::ev_to_j, 
		//	te * Constants::ev_to_j, imp_v, imp.get_mass(), me, 
		//	-imp.get_charge() * Constants::charge_e, Constants::charge_e)};
		double nu_ze {calc_momentum_loss_freq(ne, te, ne, te, imp_v, 
			imp.get_mass(), me, -imp.get_charge() * Constants::charge_e, 
			Constants::charge_e)};
		std::cout << std::scientific << "nu_ze = " << nu_ze << '\n';

		// Load ion mass in kg. Calculate momentum loss frequency due to 
		// impurity-ion collisions. Assuming singly charged ions, but the 
		// function is generalized to any charge. Just need to change what's
		// passed in for q2 here if you want other charges. Also assuming
		// quasi-neutrality and passing in ni=ne here. 
		double mi {Input::get_opt_dbl(Input::gkyl_ion_mass_amu) 
			* Constants::amu_to_kg};
		//double nu_zi {calc_momentum_loss_freq(ne, te * Constants::ev_to_j, 
		//	ti * Constants::ev_to_j, imp_v, imp.get_mass(),	mi, 
		//	-imp.get_charge() * Constants::charge_e, -Constants::charge_e)};
		double nu_zi {calc_momentum_loss_freq(ne, te, ne, ti, imp_v, 
			imp.get_mass(),	mi, -imp.get_charge() * Constants::charge_e, 
			-Constants::charge_e)};
		std::cout << std::scientific << "nu_zi = " << nu_zi << '\n';

		// Test that the equation gives expected values for e-i (just to pick
		// and example collision type). 
		double e_v {std::sqrt(3.0 * te * Constants::ev_to_j / me)};
		double nu_ei_test {1.33e5 * ne * 1.0e-20 / 
			std::pow(ti * 1.0e-3, 1.5)};
		//double nu_ei {calc_momentum_loss_freq(ne, te * Constants::ev_to_j, 
		//	ti * Constants::ev_to_j, e_v, me,	mi, 
		//	Constants::charge_e, -Constants::charge_e)};
		double nu_ei {calc_momentum_loss_freq(ne, te, ne, ti, e_v, me, mi, 
			Constants::charge_e, -Constants::charge_e)};
		std::cout << "e_v = " << e_v << '\n';
		std::cout << std::scientific << "nu_ei = " << nu_ei << " (" 
			<< nu_ei_test << ")\n";

		// Calculate total momentum (i.e., change in velocity since mass is 
		// constant) loss from each type of collision.
		// d/dt(mz*vz) = -nu * mz * vz  --> dvz = -nu * vz * dt
		double imp_dv {(nu_ze + nu_zi) * imp_v * imp_time_step};
		std::cout << "imp_dv = " << imp_dv << '\n';

		// Straightforward approach to apply this to the impurity is to scale 
		// down each impurity velocity component proportional to the amount
		// that was lost in the collisions.
		// Before: v0 = (vx0, vy0, vz0)
		// After:  v1 = a * v0 = (a*vx0, a*vy0, a*vz0) = (vx1, vy1, vz1)
		// Thus each component of the after-collision vector is:
		// vx1 = a * vx0    and likewise for y and z.
		double mom_loss_frac {(imp_v - imp_dv) / imp_v};
		std::cout << "mom_loss_frac = " << mom_loss_frac << '\n';
		std::cout << "Before: " << imp.get_vx() << ", " << imp.get_vy() << ", " 
			<< imp.get_vz() << '\n';
		imp.set_vx(imp.get_vx() * mom_loss_frac); 
		imp.set_vy(imp.get_vy() * mom_loss_frac); 
		imp.set_vz(imp.get_vz() * mom_loss_frac); 

		std::cout << "After:  " << imp.get_vx() << ", " << imp.get_vy() << ", " 
			<< imp.get_vz() << '\n';

	}

	/*
	void collision_step(Impurity::Impurity& imp, 
		const double te, const double ti, const double ne, 
		double imp_time_step)
	{

		//std::cout << "=============================\n";
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
		//std::cout << "lambda_debye = " << lambda_debye << '\n';

		// Sample a large number of different electron and ion velocites from
		// a Maxwellian distribution centered on Te, Ti. The algoritm returns
		// two numbers at a time.
		int imp_coll_num {Input::get_opt_int(Input::imp_coll_num)};

		// Calculate average ion velocity. We don't really have the standard 
		// deviation of the the ion velocity, so we're hardcoding it as avg/4.
		// This could certainly be improved in the future. Multiplying by
		// ev_to_j converts eV --> Joules.
		double avg_ion_vel {std::sqrt(2.0 * ti * Constants::ev_to_j / mi)};
		double std_ion_vel {avg_ion_vel / 4.0};
		std::vector<double> ion_vels {generate_random_vels(avg_ion_vel, 
			std_ion_vel, imp_coll_num)};
		//std::cout << "avg_ion_vel = " << avg_ion_vel << '\n';

		// Likewise for electron velocity.
		double avg_elec_vel {std::sqrt(2.0 * te * Constants::ev_to_j / me)};
		double std_elec_vel {avg_elec_vel / 4.0};
		std::vector<double> elec_vels {generate_random_vels(avg_elec_vel, 
			std_elec_vel, imp_coll_num)};
		//std::cout << "avg_elec_vel = " << avg_elec_vel << '\n';

		// Then calculate the corresponding center of mass velocities. This
		// needs the impurity velocity. We are reusing ion_vels and elec_vels
		// to save some time.
		double vz {calc_imp_vel(imp)};
		//double vz {2 * te * Constants::ev_to_j / imp.get_mass()};
		//std::cout << "vz = " << vz << '\n';
		//for (std::size_t i {}; i < ion_vels.size(); ++i)
		//{
		//	ion_vels[i] = std::sqrt(vz * vz + ion_vels[i] * ion_vels[i]);
		//}
		//for (std::size_t i {}; i < elec_vels.size(); ++i)
		//{
		//	elec_vels[i] = std::sqrt(vz * vz + elec_vels[i] * elec_vels[i]);
		//}
		std::vector<double> ion_com_vels {generate_com_vels(imp, ion_vels)};
		std::vector<double> elec_com_vels {generate_com_vels(imp, elec_vels)};

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

		// Deflection frequency for elec-imp collisions (radians/s)
		double elec_imp_freq {calc_defl_freq(red_mass_ze, elec_com_vels, elec_b, 
			Constants::charge_e, -imp.get_charge() * Constants::charge_e,
			lambda_debye, vz, imp_time_step, te, ne)};
		//std::cout << "elec_imp_ang = " << elec_imp_ang << '\n';

		// Deflection frequency for ion-imp collisions (radians/s). Assuming 
		// singly charged ions here, though the function is general and can 
		// handle other charge states if that's needed.
		double ion_imp_freq {calc_defl_freq(red_mass_zi, ion_com_vels, ion_b, 
			-Constants::charge_e, -imp.get_charge() * Constants::charge_e,
			lambda_debye, vz, imp_time_step, te, ne)};
		//std::cout << "ion_imp_ang = " << ion_imp_ang << '\n';

		// The total deflection angle is the sum of these frequencies times
		// the impurity time step.
		double tot_defl_ang {(elec_imp_freq + ion_imp_freq) * imp_time_step};

		// See how we compare with theory.
		compare_theory_90_zi(ion_imp_freq, ne, te * Constants::ev_to_j, 
			ti * Constants::ev_to_j, vz, red_mass_zi, 
			imp.get_mass(), mi, imp.get_charge() * Constants::ev_to_j, 
			Constants::ev_to_j);
		
		// At this point we have a deflection angle for an arbitrary 
		// imp_coll_num number of collisions. Over this timestep, the impurity
		// likely encounters many more collisions though. So to estimate it
		// without setting imp_coll_num to crazy large numbers, calculate by
		// what factor the real number of collisions is larger by, and multiply
		// the angle we have by that. 
	
		// Calculate how many electron/ion pairs the impurity "sees" in the
		// Debye cylinder corresponding to the impurity's distance traveled
		// during the time step (dist = vz * imp_time_step). 
		//double n_debye_cyl {Constants::pi * Constants::eps0 * te 
		//	* Constants::ev_to_j * imp_time_step * vz / (Constants::charge_e 
		//	* Constants::charge_e)};
		//std::cout << "n_debye_cyl = " << n_debye_cyl << '\n';

		//double tot_defl_ang {(n_debye_cyl / static_cast<double>(imp_coll_num)) 
		//	* (elec_imp_ang + ion_imp_ang)};	
		//std::cout << "tot_defl_ang = " << tot_defl_ang << '\n';

		// Now modify the impurity ion's trajectory by tot_defl_ang.
		//std::cout << "imp velocity\n";
		//std::cout << imp.get_vx() << ", " << imp.get_vy() << ", " 
		//	<< imp.get_vz() << '\n';
		rotate_imp(imp, tot_defl_ang);
		//std::cout << imp.get_vx() << ", " << imp.get_vy() << ", " 
		//	<< imp.get_vz() << '\n';

	}
	*/
}
