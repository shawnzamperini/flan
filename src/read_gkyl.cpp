#include <iostream>
#include <string>
/*
 * Opting out of ADIOS2 since I think it's going to get removed...
	#include <adios2.h>
*/
#include "read_input.h"

namespace Gkyl
{
	// Placeholder function
	void read_gkyl()
	{
		std::cout << "We are ready to load some Gkeyll\n";
		std::cout << "Gkeyll case name: " << Input::get_opt_str(Input::gkyl_casename) << '\n';
		std::cout << "Start frame: " << Input::get_opt_int(Input::gkyl_frame_start) << '\n';
		std::cout << "End frame: " << Input::get_opt_int(Input::gkyl_frame_end) << '\n';
	}
}
