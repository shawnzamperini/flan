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
	* @brief Update impurity velocity based on the Lorentz force
	*
	* @param imp Impurity object that is updated within function
	* @param bkg Reference to the loaded Background object
	* @param imp_time_step Size of time step specified in input file
	* @param tidx Index in time array imp is at
	* @param xidx Index in x array imp is at
	* @param yidx Index in y array imp is at
	* @param zidx Index in z array imp is at
	*/
	void lorentz_update(Impurity& imp, const Background::Background& bkg,
		const double imp_time_step, const int tidx, const int xidx, 
		const int yidx, const int zidx);
		
	/**
	* @brief Move particle based on its current velocity and the time step
	*
	* @param imp Impurity object that is updated within function
	* @param imp_time_step Size of time step specified in input file
	*/
	void step(Impurity& imp, const double imp_time_step);

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
		const Background::Background& bkg, const int tidx, const int xidx, 
		const int yidx, const int zidx);

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
	* @brief Handle impurity ion ionization and recombination
	*
	* The impurity charge is modified within this function, if necessary. If
	* the probabilities for ionization or recombination are greater than 1.0,
	* ioniz_warnings and recomb_warnings are incremented by 1. 
	*
	* @param ioniz_warnings Integer that tracks number of ionization 
	* probability > 1.0 events
	* @param recomb_warnings Integer that tracks number of recombination 
	* probability > 1.0 events
	*/
	void ioniz_recomb(Impurity& imp, const Background::Background& bkg,
		const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADAS& oa_recomb, const double imp_time_step, 
		const int tidx, const int xidx, const int yidx, const int zidx,
		int& ioniz_warnings, int& recomb_warnings);

	/**
	* @brief Controlling function for following an impurity until a terminating
	* condition is met.
	*
	* @param imp Reference to the Impurity object to follow
	* @param bkg Reference to the Background to follow the impurity in
	* @param imp_stats Reference to the Statistics object that keeps track of
	* the Monte Carlo statistics during the simulation.
	* @param oa_ioniz An OpenADAS object containing the ionization rate 
	* coefficients.
	* @param oa_recomb An OpenADAS object containing the recombination rate 
	* coefficients.
	* @param ioniz_warnings Integer that tracks number of ionization 
	* probability > 1.0 events
	* @param recomb_warnings Integer that tracks number of recombination 
	* probability > 1.0 events
	* @param imp_coll_on Boolean controlling if collisions are on or not
	*/
	void follow_impurity(Impurity& imp, const Background::Background& bkg, 
		Statistics& imp_stats, const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADAS& oa_recomb, int& ioniz_warnings, 
		int& recomb_warnings, const bool imp_coll_on);

	/**
	* @brief Print out warnings if ionization/recombination probabilities were
	* greater than 1.0 at some point.
	* @param ioniz_warnings Integer that tracks number of ionization 
	* probability > 1.0 events
	* @param recomb_warnings Integer that tracks number of recombination 
	* probability > 1.0 events
	*/
	void print_ioniz_recomb_warn(int ioniz_warnings, int recomb_warnings);
	
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
	* @param oa_ioniz An OpenADAS object containing the ionization rate 
	* coefficients.
	* @param oa_recomb An OpenADAS object containing the recombination rate 
	* coefficients.
	*/
	void main_loop(const Background::Background& bkg, Statistics& imp_stats,
		const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADAS& oa_recomb);

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
