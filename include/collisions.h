/**
* @file collisions.h
* @brief Header file for collisions.cpp
*/

#include "impurity.h"

namespace Collisions
{
	/**
	* @brief Entry point for collision model
	*/
	void collision_step(Impurity::Impurity& imp, 
		const double te, const double ti, const double ne, 
		double imp_time_step);

}
