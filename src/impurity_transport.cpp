/**
* @file  impurity_transport.cpp
*
* @brief Impurity transport routines
*/

#include <omp.h>
#include <tuple>
#include <vector>
#include <algorithm>
#include <cmath>

#include "impurity_transport.h"
#include "background.h"
#include "read_input.h"
#include "vectors.h"
#include "impurity.h"
#include "impurity_stats.h"
#include "random.h"
#include "constants.h"

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
	double get_birth_t(const Background::Background& bkg)
	{
		return Random::get(bkg.get_t_min(), bkg.get_t_max());
	}

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
	double get_birth_x(const Background::Background& bkg)
	{	
		// Load input options as local variables for cleaner code
		double xmin {Input::get_opt_dbl(Input::imp_xmin)};		
		double xmax {Input::get_opt_dbl(Input::imp_xmax)};		
	
		// Uniformily distributed between xmin, xmax.
		return Random::get(xmin, xmax);
	}

	/**
	* @brief Get starting y location to pass into an Impurity object
	*
	* Function to return a random starting x location for an impurity ion.
	* Right now this uniformily chooses between the y bounds of the simulation
	* volume.
	*
	* @param bkg Reference to the loaded Background object
	*/
	double get_birth_y(const Background::Background& bkg)
	{
		return Random::get(bkg.get_y_min(), bkg.get_y_max());
	}

	/**
	* @brief Get starting z location to pass into an Impurity object
	*
	* Function to return the starting z location for an impurity ion. Right now
	* this just uses a single value passed in via the input file.
	*
	* @param bkg Reference to the loaded Background object
	*/
	double get_birth_z(const Background::Background& bkg)
	{
		double imp_zstart_val {Input::get_opt_dbl(Input::imp_zstart_val)};
		return imp_zstart_val;
	}
	
	/**
	* @brief Get starting charge to pass into an Impurity object
	*
	* Function to return the starting charge for an impurity ion. Right now
	* this just uses a single value passed in via the input file.
	*
	* @param bkg Reference to the loaded Background object
	*/
	int get_birth_charge()
	{
		int imp_charge {Input::get_opt_int(Input::imp_init_charge)};
		return imp_charge;
	}

	/**
	* @brief Create a primary impurity ion
	*
	* Primary impurities are those that start according to the initial 
	* condition options specified in the input file. Later, these can create 
	* secondary Impurity objects as a result of Monte Carlo splitting.
	*
	* @param bkg Reference to the loaded Background object
	*/
	Impurity create_primary_imp(const Background::Background& bkg)
	{
		// Load input options as local variables for cleaner code.
		double imp_mass_amu {Input::get_opt_dbl(Input::imp_mass_amu)};

		// Starting t,x,y,z for the impurity
		double t_imp = get_birth_t(bkg);
		double x_imp = get_birth_x(bkg);
		double y_imp = get_birth_y(bkg);
		double z_imp = get_birth_z(bkg);

		// Assume starting at rest, but if this ever changes do it here
		double vx_imp {0.0};
		double vy_imp {0.0};
		double vz_imp {0.0};

		// Impurity starting charge
		int charge_imp = get_birth_charge();
	
		// Return a temporary Impurity object. Since these are primary Impurity
		// objects they by definition start with weight = 1.0.
		double imp_weight {1.0};
		double imp_mass {imp_mass_amu * Constants::amu_to_kg};
		return {t_imp, x_imp, y_imp, z_imp, vx_imp, vy_imp, vz_imp, 
			imp_weight, charge_imp, imp_mass};
	}

	/**
	* @brief Helper function to get and remove the last element in a vector
	*
	* @param vec Reference to vector to operate on
	*/
	template <typename T>
	T pop_back_remove(std::vector<T>& vec)
	{
		T element {vec.back()};
		vec.pop_back();
		return element;
	}

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
	int get_nearest_index(const std::vector<T>& vec, const T value)
	{
		auto lower = std::lower_bound(vec.begin(), vec.end(), value);
    
		if (lower == vec.begin()) 
		{
			return 0;
		}
		if (lower == vec.end()) 
		{
			return vec.size() - 1;
		}
		
		auto prev = lower - 1;
		if (std::fabs(*lower - value) < std::fabs(*prev - value)) 
		{
			return lower - vec.begin();
		} 
		else 
		{
			return prev - vec.begin();
		}
	}

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
	int get_nearest_cell_index(const std::vector<T>& grid_edges, const T value)
	{
		// Get the index of the first value in grid_edges that is larger
		// than value.
		auto lower = std::lower_bound(grid_edges.begin(), grid_edges.end(), 
			value);

		// Realize that one minus the index represented by lower is the value
		// we're after in the vectors with values at the cell centers.
		//  ____________
		//  |_0_|_1_|_2_|  <-- cell center indices
		//  0   1   2   3  <-- grid_edges indices
		//          ^
		//        lower
		//
		// In this example, we want 1 returned, so we return 2 - 1 = 1. 
		if (lower == grid_edges.begin()) return 0;
		else return lower - grid_edges.begin() - 1;
	}

	void do_lorentz_step(Impurity& imp, const Background::Background& bkg,
		const double imp_time_step, const int tidx, const int xidx, 
		const int yidx, const int zidx)
	{
		// Impurity's charge
		double imp_q {imp.get_charge() * Constants::charge_e};

		//std::cout << "Ex, Ey, Ez = " << bkg.get_ex()(tidx, xidx, yidx, zidx)
		//	<< bkg.get_ey()(tidx, xidx, yidx, zidx) 
		//	<< bkg.get_ez()(tidx, xidx, yidx, zidx) << "\n";
		//std::cout << "Bz = " << bkg.get_b()(tidx, xidx, yidx, zidx) << "\n";

		// Each component of the Lorentz force
		double fx {imp_q * (bkg.get_ex()(tidx, xidx, yidx, zidx) + 
			imp.get_vy() * bkg.get_b()(tidx, xidx, yidx, zidx))};
		double fy {imp_q * (bkg.get_ey()(tidx, xidx, yidx, zidx) + 
			-imp.get_vx() * bkg.get_b()(tidx, xidx, yidx, zidx))};
		double fz {imp_q * bkg.get_ez()(tidx, xidx, yidx, zidx)};

		// Change in velocity over time step
		double dvx {fx * imp_time_step / imp.get_mass()};
		double dvy {fy * imp_time_step / imp.get_mass()};
		double dvz {fz * imp_time_step / imp.get_mass()};

		// Update particle velocity
		imp.set_vx(imp.get_vx() + dvx);
		imp.set_vy(imp.get_vy() + dvy);
		imp.set_vz(imp.get_vz() + dvz);

		// Update particle position
		imp.set_x(imp.get_x() + imp.get_vx() * imp_time_step);
		imp.set_y(imp.get_y() + imp.get_vy() * imp_time_step);
		imp.set_z(imp.get_z() + imp.get_vz() * imp_time_step);
	}

	void record_stats(Statistics& imp_stats, const Impurity& imp, const int tidx, 
		const int xidx, const int yidx, const int zidx)
	{
		// Add one to counts to this location	
		imp_stats.add_counts(tidx, xidx, yidx, zidx, 1);

		// Add particle's weight to this location. Note: We should be adding
		// the particle weight * imp_time_step, but since the time step is
		// constant, we can add it at the end and save a floating point
		// operation.
		imp_stats.add_weights(tidx, xidx, yidx, zidx, imp.get_weight());
	}

	bool check_boundary(const Background::Background& bkg, Impurity& imp)
	{
		// The x boundaries are treated as absorbing. Could certianly add
		// different options here in the future.
		if (imp.get_x() <= bkg.get_grid_x()[0])
		{
			//std::cout << "absorbed: lower_x " << imp.get_x() << " < " 
			//	<< bkg.get_grid_x()[0] << "\n";
			return false;
		}
		if (imp.get_x() >= bkg.get_grid_x().back()) 
		{
			//std::cout << "absorbed: upper_x " << imp.get_x() << " > " 
			//	<< bkg.get_grid_x().back() << "\n";
			return false;
		}

		// Absorbing z boundaries
		if (imp.get_z() <= bkg.get_grid_z()[0]) return false;
		if (imp.get_z() >= bkg.get_grid_z().back()) return false;

		// Periodic y boundaries. This means if an impurity crosses one of
		// these boundaries it will wrap around to the other one.
		if (imp.get_y() < bkg.get_grid_y()[0])
		{
			double dy {imp.get_y() - bkg.get_grid_y()[0]};
			imp.set_y(bkg.get_grid_y().back() + dy);
		}
		if (imp.get_y() > bkg.get_grid_y().back())
		{
			double dy {imp.get_y() - bkg.get_grid_y().back()};
			imp.set_y(bkg.get_grid_y()[0] + dy);
		}

		return true;
	}

	void follow_impurity(Impurity& imp, const Background::Background& bkg, 
		Statistics& imp_stats)
	{
		// Timestep of impurity transport simulation
		double imp_time_step {Input::get_opt_dbl(Input::imp_time_step)};

		bool continue_following {true};
		while (continue_following)
		{
			// Debugging
			//std::cout << "imp t, x, y, z: " << imp.get_t() << ", " << 
			//	imp.get_x() << ", " << imp.get_y() << ", " << imp.get_z() 
			//	<< '\n';

			// Get nearest time index
			int tidx {get_nearest_index(bkg.get_times(), imp.get_t())};

			// Get x,y,z indices
			int xidx {get_nearest_cell_index(bkg.get_grid_x(), imp.get_x())};
			int yidx {get_nearest_cell_index(bkg.get_grid_y(), imp.get_y())};
			int zidx {get_nearest_cell_index(bkg.get_grid_z(), imp.get_z())};

			// Perform a step according to the Lorentz force
			do_lorentz_step(imp, bkg, imp_time_step, tidx, xidx, yidx, zidx);

			// Update particle time
			imp.set_t(imp.get_t() + imp_time_step);

			// Update statistics
			record_stats(imp_stats, imp, tidx, xidx, yidx, zidx);

			// Check for a terminating or boundary conditions
			continue_following = check_boundary(bkg, imp);

			// Check for a collision

			// Check for ionization or recombination
		}

	}

	// Main particle following loop. This is where the magic begins.
	void main_loop(const Background::Background& bkg, Statistics& imp_stats)
	{
		// Load input options as local variables for cleaner code.
		int imp_num {Input::get_opt_int(Input::imp_num)};

		// https://stackoverflow.com/questions/29633531/user-defined-
		// reduction-on-vector-of-varying-size/29660244#29660244
		#pragma omp declare reduction(+: Statistics: \
			omp_out = omp_out + omp_in) \
			initializer(omp_priv(omp_orig))

		std::cout << "Starting particle following...\n";
	
		// Loop through one impurity at a time, tracking it from its birth
		// time/location to the end. 
		// dynamic scheduling likely the best here since the loop times can
		// vary widely.
		// Probably want to make this a while loop with OMP tasks for logic
		// that goes like this (serial example):
		//   split_imps {}
		//   do
		//     *follow impurity routine*
		//     *potentially add imp to split_imps*
		//   while (!split_imps.is_empty())
		int thread_imp_count {};
		int thread_imp_num {};
		#pragma omp parallel for schedule(dynamic) \
			shared(bkg) \
			firstprivate(thread_imp_num, thread_imp_count) \
			reduction(+: imp_stats)
		for (int i = 0; i < imp_num; ++i)
		{
			// OpenMP overhead
			[[maybe_unused]] int thread_id {omp_get_thread_num()};
			[[maybe_unused]] int num_threads {omp_get_num_threads()};
			thread_imp_num = imp_num / num_threads;

			// Printout of progress. We don't want to use atomic or anything
			// that can slow the program down, so we will just report progress
			// for a single thread since, in theory, all threads should finish
			// near the same time. For smaller numbers of particles this can
			// give silly numbers, but it should give the correct numbers at
			// production level numbers of impurities.
			if (thread_id == 0)
			{
				// Print progress this many times
				int prog_interval {10};

				// Percent that we've followed
				double perc_complete {static_cast<double>(thread_imp_count) 
					/ thread_imp_num * 100};

				// Only do this is enough impurities are being followed, 
				// otherwise the following modulo in the if statement will
				// be a divide by zero error.
				if (thread_imp_num > prog_interval)
				{
					if ((thread_imp_count % (thread_imp_num / prog_interval)) == 0 
						&& thread_imp_count > 0)
					{
						std::cout << "Followed " << thread_imp_count * num_threads 
							<< "/" << imp_num << " impurities (" 
							<< static_cast<int>(perc_complete) << "%)\n";
					}
				}
			}

			// Create starting impurity ion
			Impurity primary_imp = create_primary_imp(bkg);
			
			std::vector<Impurity> imps {primary_imp};
			while (!imps.empty())
			{
				Impurity imp {pop_back_remove(imps)};

				follow_impurity(imp, bkg, imp_stats);
			}

			// Increment thread-specific counter. This just counts primary
			// impurities. 
			thread_imp_count += 1;
		}
	}

	// Entry point to particle following routine(s).
	Statistics follow_impurities(Background::Background& bkg)
	{
		// Initialize particle statistics vectors, all contained within a
		// Statistics object.
		Statistics imp_stats {bkg.get_dim1(), bkg.get_dim2(), 
			bkg.get_dim3(), bkg.get_dim4()};

		// Execute main particle following loop.
		main_loop(bkg, imp_stats);
	
		// Convert the statistics into meaningful quantities. We are scaling
		// the density by the scale factor * time step. Technically speaking,
		// the time step should be accounted for in the main loop each time
		// we score a particle (weight * time step instead of just weight), 
		// but since the time step is constant, we can avoid that floating 
		// point operation in the main loop and apply it here after the fact.
		std::cout << "Calculating derived quantities...\n";
		int imp_num {Input::get_opt_int(Input::imp_num)};
		double imp_source_scale_fact {
			Input::get_opt_dbl(Input::imp_source_scale_fact)};
		double imp_time_step {Input::get_opt_dbl(Input::imp_time_step)};
		imp_stats.calc_density(bkg, imp_num, 
			imp_source_scale_fact * imp_time_step);
		
		return imp_stats;

	}

}
