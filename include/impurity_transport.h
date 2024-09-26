/**
* @file impurity_transport.h
*
* @brief Header file for impurity_transport.cpp
*/

#ifndef IMPURITY_TRANSPORT_H
#define IMPURITY_TRANSPORT_H

#include "background.h"
#include "impurity.h"
#include "impurity_stats.h"

namespace Impurity
{
	// Functions related to intializing an Impurity object.
	double get_birth_t(const Background::Background& bkg);
	double get_birth_x(const Background::Background& bkg);
	double get_birth_y(const Background::Background& bkg);
	double get_birth_z(const Background::Background& bkg);
	int get_birth_charge();
	Impurity create_primary_imp(const Background::Background& bkg);
	
	// Helper function to return an element from a vector while also removing
	// it.
	template <typename T>
	T pop_back_remove(std::vector<T>& vec);

	// Returns the nearest index in a vector to a value.
	template <typename T>
	int get_nearest_index(const std::vector<T>& vec, const T value);

	// Return the nearest index in the cell center arrays (e.g., those in the 
	// Background object) to value. This is quicker than get_nearest_index
	// because it does not determine which edge it is closer to. grid_edges
	// must be sorted.
	template <typename T>
	int get_nearest_cell_index(const std::vector<T>& grid_edges, const T value);

	// Perform an impurity step based on the Lorentz force.
	void do_lorentz_step(Impurity& imp, const Background::Background& bkg,
		const double imp_time_step, const int tidx, const int xidx, 
		const int yidx, const int zidx);
		
	// Record Monte Carlo statistics for imp into imp_stats
	void record_stats(Statistics& imp_stats, const Impurity& imp, 
		const int tidx, const int xidx, const int yidx, const int zidx);

	// Check if an Impurity has encountered a boundary condition. Returns
	// false if following is to be finished, true if otherwise.
	bool check_boundary(const Background::Background& bkg, Impurity& imp);
	
	// Follow a single Impurity until a terminating condition is met.
	void follow_impurity(Impurity& imp, const Background::Background& bkg, 
		Statistics& imp_stats);

	// Main loop for impurity transport simulation.
	void main_loop(const Background::Background& bkg, Statistics& imp_stats);

	// Entry level function for impurity transport simulation.
	Statistics follow_impurities(Background::Background& bkg);
}

#endif
