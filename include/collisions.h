/**
* @file collisions.h
* @brief Header file for collisions.cpp
*/

#include "impurity.h"
#include "options.h"
#include "background.h"

namespace Collisions
{
	/**
	* @brief Update impurity Cartesian velocity according to Nanbu collision
	* model.
	*/
	void nanbu_coll(Impurity::Impurity& imp, const Background::Background& bkg,
		const int tidx, const int xidx, const int yidx, const int zidx,
		const Options::Options& opts, bool elec, const double imp_time_step);

}
