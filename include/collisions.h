/**
* @file collisions.h
* @brief Header file for collisions.cpp
*/

#include "impurity.h"
#include "background.h"

namespace Collisions
{
	/**
	* @brief Calculate variable time step for reasonable collision calculation
	*/
	void set_var_time_step_coll(double& dt_coll, const Impurity::Impurity& imp,
		const Background::Background& bkg, const int tidx, const int xidx, 
		const int yidx, const int zidx);
		
	/**
	* @brief Entry point for collision model
	*/
	void collision_update(Impurity::Impurity& imp, 
		const double te, const double ti, const double ne, 
		double imp_time_step);

}
