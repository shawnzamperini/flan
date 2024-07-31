#include <iostream>
#include <fstream>
#include <string>

#include "read_input.h"
#include "read_gkyl.h"
#include "impurity_transport.h"

int main(int argc, char *argv[])
{
	std::cout << "Welcome to Flan V X.X\n";

	// Read in input options, overwriting all the default values.
	if (argc != 2)
	{
		std::cerr << "Error! Enter case name after flan. Example:\n";
		std::cerr << "  ./flan case_name\n";
		return -1;
	}

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
	Gkyl::Background bkg {Gkyl::read_gkyl()};
	
	// Interpolate additional frames between each Gkeyll frame to artificially
	// increase the time resolution of the simulation.
	// To-do: Implement this later once I have an easier way to test it.
	

	// Begin main particle following loop.
	Impurity::follow_impurities(bkg);

	// Save simulation results.
	

	return 0;
}
