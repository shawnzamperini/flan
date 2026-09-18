/**
* @file collisions.h
* @brief Header file for collisions.cpp
*/

#include "background.h"
#include "impurity.h"
#include "impurity_stats.h"
#include "options.h"
#include "pcg32.h"
#include "slots.h"

namespace Collisions
{
	/**
	* @brief Update impurity Cartesian velocity according to Nanbu collision
	* model.
	*/
	void nanbu_coll(Slots::Slots& slots, const Background::Background& bkg,
		const Options::Options& opts, bool elec, const double dt, 
		Impurity::Statistics& imp_stats, std::vector<pcg32>& rngs);
}
