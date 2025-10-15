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
