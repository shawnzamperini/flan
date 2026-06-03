#include "background.h"
#include "mpi.h"
#include "options.h"
#include "read_gkyl.h"
#include "read_test.h"

namespace Background
{

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

		// Unrecognized option
		else
		{
			std::cerr << "Error! bkg_source = " << opts.bkg_source() << " not"
				<< " recognized.\n";
		}

		return bkg;
	}

	
}
