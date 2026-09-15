/**
* @file collisions.h
* @brief Header file for collisions.cpp
*/

#include "impurity.h"
#include "impurity_stats.h"
#include "options.h"
#include "background.h"
#include "slots.h"

namespace Collisions
{
	/**
	* @brief Update impurity Cartesian velocity according to Nanbu collision
	* model.
	*/
	void nanbu_coll(Slots::Slots& slots, const Background::Background& bkg,
		const Options::Options& opts, bool elec, const double dt);
}
