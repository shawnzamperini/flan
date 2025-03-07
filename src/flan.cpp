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

void flan(const Inputs& inpts)
{
	// Load version info from config.h file
	std::cout << "Welcome to Flan v" << PROJECT_VERSION << "\n";

	// Load options
	Options::Options opts {Options::load_options(inpts)};

	// Load in a Gkeyll run for the background plasma, returns a Background
	// class object.
	Background::Background bkg {Gkyl::read_gkyl(opts)};

	// Interpolate additional frames between each Gkeyll frame to artificially
	// increase the time resolution of the simulation.
	// To-do: Implement this later once I have an easier way to test it.
	
	// Begin main particle following loop.
	Impurity::Statistics imp_stats {Impurity::follow_impurities(bkg, opts)};

	// Save simulation results.
	SaveResults::save_results(bkg, imp_stats, opts);

}
