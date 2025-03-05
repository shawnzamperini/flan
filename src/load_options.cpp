#include <iostream>

#include "load_options.h"
#include "options.h"

namespace Options
{
	Options load_options()
	{
		// Object to hold all the options
		Options opts {};
		
		// First load in all the defaults, then overwrite the values with
		// whatever the user input
		load_defaults();
		load_input();

		return opts;
	}

	void load_defaults()
	{
		std::cout << "Loading defaults...\n";

	}

	void load_input()
	{
		std::cout << "Loading input...\n";

	}

}
