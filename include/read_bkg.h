#include "background.h"
#include "options.h"

namespace Background
{

	/**
	* @brief Read in plasma background
	*
	* @param opts The Options object controlling the simulation
	*
	* @returns A Background object
	*/
	Background read_bkg(const Options::Options& opts);
}
