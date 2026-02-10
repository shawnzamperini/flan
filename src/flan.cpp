/**
* @file  flan.cpp
*
* @brief Entry point and top level control of Flan library.
*/

#include <iostream>
#include <fstream>
#include <string>

#include "background.h"
#include "config.h"
#include "flan.h"
#include "flan_types.h"
#include "impurity_stats.h"
#include "impurity_transport.h"
#include "load_options.h"
#include "options.h"
#include "read_gkyl.h"
#include "save_results.h"
#include "timer.h"
#include "vectors.h"

void flan(const Inputs& inpts)
{
	// Create Timer object and begin total simulation timer
	Timer::Timer timer {};
	timer.start_total_timer();

	// For printing output, flush buffer after every output operation
	std::cout << std::unitbuf;

	// Load version info from config.h file
	std::cout << "Welcome to Flan v" << PROJECT_VERSION << '\n';
	std::cout << "Commit: " << GIT_COMMIT_HASH << '\n';

	// Load options
	Options::Options opts {Options::load_options(inpts)};

	// Load in a Gkeyll run for the background plasma, returns a Background
	// class object.
	timer.start_read_timer();
	Background::Background bkg {Gkyl::read_gkyl(opts)};
	timer.end_read_timer();

	/*
	// DEBUG: Flan is missing a test suite. For the paper, the reviewer wants
	// some kind of test, so I am going to overwrite the background with 
	// constant values to test simple particle gyromotion.
	std::cout << "WARNING! Replacing background with constant values for debugging!\n";

	// Vector to hold constant values that have X,Y,Z components
	std::vector<BkgFPType> const_eX_t0 (std::ssize(bkg.get_z()));
	std::vector<BkgFPType> const_eY_t0 (std::ssize(bkg.get_z()));
	std::vector<BkgFPType> const_eZ_t0 (std::ssize(bkg.get_z()));
	std::vector<BkgFPType> const_bX_t0 (std::ssize(bkg.get_z()));
	std::vector<BkgFPType> const_bY_t0 (std::ssize(bkg.get_z()));
	std::vector<BkgFPType> const_bZ_t0 (std::ssize(bkg.get_z()));
	for (int l {}; l < std::ssize(bkg.get_z()); ++l)
	{
		// Values at the middle of x, y, z that we will set for the constant 
		// values. Using midpoint in x, y directions.
		const_eX_t0[l] = 0.0; 
		const_eY_t0[l] = 0.0; 
		if (false)
		{
			const_eZ_t0[l] = 0.0; 
		}
		else
		{
			const_eZ_t0[l] = bkg.get_eZ()(0, 48, 32, l); 
		}
		const_bX_t0[l] = bkg.get_bX()(0, 48, 32, l); 
		const_bY_t0[l] = bkg.get_bY()(0, 48, 32, l); 
		const_bZ_t0[l] = bkg.get_bZ()(0, 48, 32, l); 

		// Useful printout
		std::cout << "---------------------------\n";
		std::cout << "zidx = " << l << '\n';
		std::cout << "  EZ = " << const_eZ_t0[l] << '\n';
		std::cout << "  BX = " << const_bX_t0[l] << '\n';
		std::cout << "  BY = " << const_bY_t0[l] << '\n';
		std::cout << "  BZ = " << const_bZ_t0[l] << '\n';
	}

	#pragma omp parallel for
	for (int i = 0; i < std::ssize(bkg.get_times()); ++i)  // t
	{
	for (int j {}; j < std::ssize(bkg.get_x()); ++j)  // x
	{
	for (int k {}; k < std::ssize(bkg.get_y()); ++k)  // y
	{
	for (int l {}; l < std::ssize(bkg.get_z()); ++l)  // z
	{

		// Calculate index in 1D vector
		int idx {bkg.get_ne().calc_index(i,j,k,l)};

		// Replace with constant values
		bkg.get_ne().get_data()[idx] = 1e18;
		bkg.get_te().get_data()[idx] = 10;
		bkg.get_ti().get_data()[idx] = 10;
		bkg.get_vp().get_data()[idx] = 0.0; // Not needed

		// Need to account for toroidal curvature.
		bkg.get_eX().get_data()[idx] = const_eX_t0[l];
		bkg.get_eY().get_data()[idx] = const_eY_t0[l];
		bkg.get_eZ().get_data()[idx] = const_eZ_t0[l];
		bkg.get_bX().get_data()[idx] = const_bX_t0[l];
		bkg.get_bY().get_data()[idx] = const_bY_t0[l];
		bkg.get_bZ().get_data()[idx] = const_bZ_t0[l];

		// Shouldn't matter without collisions but set to zero anyways
		bkg.get_uX().get_data()[idx] = 0.0;
		bkg.get_uY().get_data()[idx] = 0.0;
		bkg.get_uZ().get_data()[idx] = 0.0;

	}
	}
	}
	}
	*/

	// Interpolate additional frames between each Gkeyll frame to artificially
	// increase the time resolution of the simulation.
	// To-do: Implement this later once I have an easier way to test it.
	
	// Begin main particle following loop.
	timer.start_imp_timer();
	Impurity::Statistics imp_stats {Impurity::follow_impurities(bkg, opts, 
		timer)};
	timer.end_imp_timer();

	// Save simulation results.
	timer.start_save_timer();
	SaveResults::save_results(bkg, imp_stats, opts);
	timer.end_save_timer();

	// End timer, print summary
	timer.end_total_timer();
	timer.print_summary();

}
