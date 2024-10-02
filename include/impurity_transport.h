/**
* @file impurity_transport.h
*
* @brief Header file for impurity_transport.cpp
*/

#ifndef IMPURITY_TRANSPORT_H
#define IMPURITY_TRANSPORT_H

#include <vector> 

#include "background.h"
#include "impurity.h"
#include "impurity_stats.h"
#include "openadas.h"


namespace Impurity
{
	/**
	* @brief Get starting time to pass into an Impurity object.
	*
	* Function to return a random starting time for an impurity ion. This is
	* obviously very basic, to have as its own function, but we do this since
	* this will probably get more complicated in the future.
	*
	* @param bkg Reference to the loaded Background object
	*/
	double get_birth_t(const Background::Background& bkg);

	/**
	* @brief Get starting x location to pass into an Impurity object
	*
	* Function to return a random starting x location for an impurity ion.
	* Right now this uniformily chooses between two values given by the user
	* in the input file. If both those values are the same it just returns
	* that number.
	*
	* @param bkg Reference to the loaded Background object
	*/
	double get_birth_x(const Background::Background& bkg);

	/**
	* @brief Get starting y location to pass into an Impurity object
	*
	* Function to return a random starting x location for an impurity ion.
	* Right now this uniformily chooses between the y bounds of the simulation
	* volume.
	*
	* @param bkg Reference to the loaded Background object
	*/
	double get_birth_y(const Background::Background& bkg);

	/**
	* @brief Get starting z location to pass into an Impurity object
	*
	* Function to return the starting z location for an impurity ion. Right now
	* this just uses a single value passed in via the input file.
	*
	* @param bkg Reference to the loaded Background object
	*/
	double get_birth_z(const Background::Background& bkg);

	/**
	* @brief Get starting charge to pass into an Impurity object
	*
	* Function to return the starting charge for an impurity ion. Right now
	* this just uses a single value passed in via the input file.
	*/
	int get_birth_charge();

	/**
	* @brief Create a primary impurity ion
	*
	* Primary impurities are those that start according to the initial 
	* condition options specified in the input file. Later, these can create 
	* secondary Impurity objects as a result of Monte Carlo splitting.
	*
	* @param bkg Reference to the loaded Background object
	*/
	Impurity create_primary_imp(const Background::Background& bkg);
	
	/**
	* @brief Helper function to get and remove the last element in a vector
	*
	* @param vec Reference to vector to operate on
	*/
	template <typename T>
	T pop_back_remove(std::vector<T>& vec);

	/**
	* @brief Function that return the nearest index in a vector to value
	*
	* This is a slightly tricky algorithm, so I relied on Microsoft's Copilot
	* AI to help me figure this one out using the standard library so it would
	* be fast.
	*
	* @param vec The vector to search in
	* @param value The value we are looking for the nearest index to
	*/
	template <typename T>
	int get_nearest_index(const std::vector<T>& vec, const T value);

	/**
	* @brief Find nearest cell index in a vector representing a grid
	*
	* This algorithm is a bit quicker than get_nearest_index because it uses
	* the grid edges. In this case, we do not care which grid edge we are
	* closer to. As long as grid_edges is sorted, once we know the index of
	* the first element that is larger than value, we have the cell index. That
	* lets us skip a few extra step that get_nearest_index has to do.
	*
	* @param grid_edges Vector containing the grid edges for the cells
	* @param value Value to find nearest cell for
	*/
	template <typename T>
	int get_nearest_cell_index(const std::vector<T>& grid_edges, const T value);

	/**
	* @brief Perform impurity step based on the Lorentz force
	*
	* @param imp Impurity object that is updated within function
	* @param bkg Reference to the loaded Background object
	* @param imp_time_step Size of time step specified in input file
	* @param tidx Index in time array imp is at
	* @param xidx Index in x array imp is at
	* @param yidx Index in y array imp is at
	* @param zidx Index in z array imp is at
	*/
	void do_lorentz_step(Impurity& imp, const Background::Background& bkg,
		const double imp_time_step, const int tidx, const int xidx, 
		const int yidx, const int zidx);
		
	/**
	* @brief Record/score particle in the ImpurityStats object
	*
	* @param imp_stats A Statistics object tracking the Monte Carlo statistics
	* of this simulation.
	* @param imp The Impurity object to score
	* @param tidx Index in time array imp is at
	* @param xidx Index in x array imp is at
	* @param yidx Index in y array imp is at
	* @param zidx Index in z array imp is at
	*/
	void record_stats(Statistics& imp_stats, const Impurity& imp, 
		const int tidx, const int xidx, const int yidx, const int zidx);

	/**
	* @brief Check if an Impurity has encountered a boundary condition.
	*
	* A boundary condition could be absorbing or reflecting, depends on the 
	* simulation settings and boundary.
	*
	* @param bkg Reference to the loaded Background object
	* @param imp The Impurity object to check for boundary condition
	*
	* @return Returns false if particle following is to be terminated, true
	* otherwise.
	*/
	bool check_boundary(const Background::Background& bkg, Impurity& imp);
	
	/**
	* @brief Controlling function for following an impurity until a terminating
	* condition is met.
	*
	* @param imp Reference to the Impurity object to follow
	* @param bkg Reference to the Background to follow the impurity in
	* @param imp_stats Reference to the Statistics object that keeps track of
	* the Monte Carlo statistics during the simulation.
	*/
	void follow_impurity(Impurity& imp, const Background::Background& bkg, 
		Statistics& imp_stats);

	/**
	* @brief Main particle following loop for impurity transport
	*
	* The is the main function for the impurity transport part of Flan.
	* The main loop is parallelized with OpenMP. Each thread is assigned an
	* impurity to follow, including any secondary impurities spawned by that
	* impurity. Dyanamic scheduling ensures each thread stays busy.
	*
	* @param bkg Reference to Background to perform simulation in
	* @param imp_stats Reference to the Statistics object that keeps track of
	* the Monte Carlo statistics during the simulation.
	*/
	void main_loop(const Background::Background& bkg, Statistics& imp_stats);

	/**
	* @brief Entry level function to impurity following routines
	*
	* This function controls the top-level functions of the impurity transport
	* simulation. See main_loop for more information.
	*
	* @param bkg Reference to Background to perform simulation in
	*
	* @return Return a Statistics object with the impurity-related results of
	* the simulation.
	*/
	Statistics follow_impurities(Background::Background& bkg);
}

#endif
