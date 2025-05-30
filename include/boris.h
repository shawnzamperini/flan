/**
* @file boris.h
*
* @brief Header file for boris.cpp.
*/

#ifndef BORIS_H
#define BORIS_H

#include "background.h"
#include "impurity.h"


namespace Boris
{
	/**
	* @brief Update impurity velocity using the Boris algorithm.
	*
	* A byproduct of this algorithm is that the particle velocity that is
	* stored in the Impurity object is actually the velocity at t - dt/2.
	*/
	void update_velocity(Impurity::Impurity& imp, 
		const Background::Background& bkg, const double dt, 
		const int tidx, const int xidx, const int yidx, const int zidx);
}

#endif 
