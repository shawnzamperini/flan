/**
* @file  flan_main.cpp
*
* @brief Entry point and main top level control of Flan.
*/

#include <iostream>
#include <fstream>
#include <string>

#include "read_input.h"
#include "read_gkyl.h"
#include "impurity_transport.h"
#include "save_results.h"
#include "background.h"
#include "impurity_stats.h"
#include "config.h"

int main(int argc, char *argv[])
{
	// Load version info from config.h file
	std::cout << "Welcome to Flan v" << PROJECT_VERSION << "\n";

	// Read in input options, overwriting all the default values.
	if (argc != 2)
	{
		std::cerr << "Error! Enter case name after flan. Example:\n";
		std::cerr << "  ./flan case_name\n";
		return -1;
	}

	// Test that we won't run into any NetCDF (really, HDF5) errors when we
	// want to save at the end. Annoying if it happens.
	// To-do.

	// Load input file.
	std::string case_name {argv[1]};
	std::ifstream input_stream {case_name + ".flin"};

	// If we couldn't open the input file, exit.
	if (!input_stream)
	{
		std::cerr << "Unable to open file. Make sure filename uses the .flin"
			" extension and call the code with just the case name. Example:\n";
		std::cerr << "  ./flan case_name\n";
		return -1;
	}
	
	// Entry function to load all the options into classes (e.g., string
	// options are in OptionStr classes, int options in OptionInt, etc.).
	Input::load_input_opts(input_stream);
	
	// Load in a Gkeyll run for the background plasma, returns a Background
	// class object.
	Background::Background bkg {Gkyl::read_gkyl()};
	
	// Interpolate additional frames between each Gkeyll frame to artificially
	// increase the time resolution of the simulation.
	// To-do: Implement this later once I have an easier way to test it.
	

	// Begin main particle following loop.
	Impurity::Statistics imp_stats {Impurity::follow_impurities(bkg)};

	// Save simulation results.
	SaveResults::save_results(case_name, bkg, imp_stats);

	return 0;
}
