/**
* @file  read_bkg.cpp
*
* @brief Entry point for loading the plasma background, which dispatches to
*        the reader for the requested background source.
*/
#include <iostream>

#include "background.h"
#include "options.h"
#include "read_gkyl.h"
#include "read_test.h"

namespace Background
{
	/**
	* Reads in the plasma background from the source chosen by the
	* bkg_source option.
	*
	* "test" loads a hardcoded test background with Test::read_test, and
	* "gkeyll" loads a Gkeyll simulation with Gkyl::read_gkyl. SOLEDGE3X
	* backgrounds (control integer 2) are not implemented yet. If the source
	* is not implemented or not recognized, an error is printed and an empty
	* Background is returned.
	*
	* @param opts Simulation options.
	* @return Background object holding the plasma background.
	*/
	Background read_bkg(const Options::Options& opts)
	{
		// Background object to fill in
		Background bkg {};

		// Load in hardcoded plasma background for testing
		if (opts.bkg_source_int() == 0)
		{
			bkg = Test::read_test(opts);
		}

		// Load in a Gkeyll background
		else if (opts.bkg_source_int() == 1)
		{
			bkg = Gkyl::read_gkyl(opts);
		}

		// Load in a SOLEDGE3X background
		else if (opts.bkg_source_int() == 2)
		{
			std::cerr << "Error! SOLEDGE3X backgrounds are not "
				<< "implemented yet!\n";
		}

		// Unrecognized option
		else
		{
			std::cerr << "Error! bkg_source = " << opts.bkg_source() << " not"
				<< " recognized.\n";
		}

		return bkg;
	}  // read_bkg
}  // namespace Background
