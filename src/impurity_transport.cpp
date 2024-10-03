/**
* @file impurity_transport.cpp
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
#include "openadas.h"


namespace Impurity
{
	double get_birth_t(const Background::Background& bkg)
	{
		return Random::get(bkg.get_t_min(), bkg.get_t_max());
	}

	double get_birth_x(const Background::Background& bkg)
	{	
		// Load input options as local variables for cleaner code
		double xmin {Input::get_opt_dbl(Input::imp_xmin)};		
		double xmax {Input::get_opt_dbl(Input::imp_xmax)};		
	
		// Uniformily distributed between xmin, xmax.
		return Random::get(xmin, xmax);
	}

	double get_birth_y(const Background::Background& bkg)
	{
		return Random::get(bkg.get_y_min(), bkg.get_y_max());
	}

	double get_birth_z(const Background::Background& bkg)
	{
		double imp_zstart_val {Input::get_opt_dbl(Input::imp_zstart_val)};
		return imp_zstart_val;
	}
	
	int get_birth_charge()
	{
		int imp_charge {Input::get_opt_int(Input::imp_init_charge)};
		return imp_charge;
	}

	Impurity create_primary_imp(const Background::Background& bkg)
	{
		// Load input options as local variables for cleaner code.
		double imp_mass_amu {Input::get_opt_dbl(Input::imp_mass_amu)};
		int imp_atom_num {Input::get_opt_int(Input::imp_atom_num)};

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
			imp_weight, charge_imp, imp_mass, imp_atom_num};
	}

	template <typename T>
	T pop_back_remove(std::vector<T>& vec)
	{
		T element {vec.back()};
		vec.pop_back();
		return element;
	}

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

	void record_stats(Statistics& imp_stats, const Impurity& imp, 
		const int tidx, const int xidx, const int yidx, const int zidx)
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
			return false;
		}
		if (imp.get_x() >= bkg.get_grid_x().back()) 
		{
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

	void ioniz_recomb(Impurity& imp, const Background::Background& bkg,
		const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADAS& oa_recomb, const double imp_time_step, 
		const int tidx, const int xidx, const int yidx, const int zidx, 
		int& ioniz_warnings, int& recomb_warnings)
	{
		// We use ne and te more than once here, so to avoid indexing
		// multiple times it is cheaper to just do it once and save it in
		// a local variable.
		double local_ne = bkg.get_ne()(tidx, xidx, yidx, zidx);
		double local_te = bkg.get_te()(tidx, xidx, yidx, zidx);
		
		// Ionization rate coefficients are indexed by charge. This is 
		// because the zeroeth charge index in the underlying rate data is
		// for neutral ionization (charge = 0), W0 --> W1+.
		double ioniz_rate {};
		if (imp.get_charge() < imp.get_atom_num())
		{
			ioniz_rate = oa_ioniz.get_rate_coeff(imp.get_charge(), 
				local_ne, local_te);
		}

		// Recombination rate coefficients are indexed by charge-1. This is
		// because this zeroeth entry is for W1+ --> W0. So if we want that
		// rate coefficient for, say, W1+, we need to pass it charge-1 so it
		// chooses that zeroeth index.
		double recomb_rate {};
		if (imp.get_charge() > 0)
		{
			recomb_rate = oa_recomb.get_rate_coeff(imp.get_charge()-1, 
			local_ne, local_te);
		}

		/*
		std::cout << "-------------------------------\n";
		std::cout << "te = " << bkg.get_te()(tidx, xidx, yidx, zidx) << '\n';
		std::cout << "ne = " << bkg.get_ne()(tidx, xidx, yidx, zidx) << '\n';
		std::cout << "charge = " << imp.get_charge() << '\n';
		std::cout << "ioniz_rate =  " << ioniz_rate << '\n';
		std::cout << "recomb_rate = " << recomb_rate << '\n';
		std::cout << "-------------------------------\n";
		*/

		// The probability of either ionization or recombination occuring is:
		//   prob = rate [m3/s] * ne [m-3] * dt [s]
		double ioniz_prob {ioniz_rate * local_ne * imp_time_step};
		double recomb_prob {recomb_rate * local_ne * imp_time_step};

		// Track number of times the probabilities are greater than 1. This
		// indicates that a smaller timestep should be used if there are a
		// significant number of warnings. 
		if (ioniz_prob > 1.0) ioniz_warnings += 1;
		if (recomb_prob > 1.0) recomb_warnings += 1;

		// For each process, pull a random number. If that number is less than
		// prob, then that event occurs. If both events occur, then they just 
		// cancel each other out and there's no change.
		if (Random::get(0.0, 1.0) < ioniz_prob)
		{
			imp.set_charge(imp.get_charge() + 1);
		}
		if (Random::get(0.0, 1.0) < recomb_prob)
		{
			imp.set_charge(imp.get_charge() - 1);
		}
	}

	void follow_impurity(Impurity& imp, const Background::Background& bkg, 
		Statistics& imp_stats, const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADAS& oa_recomb, int& ioniz_warnings, 
		int& recomb_warnings)
	{
		// Timestep of impurity transport simulation
		double imp_time_step {Input::get_opt_dbl(Input::imp_time_step)};

		bool continue_following {true};
		while (continue_following)
		{
			// Debugging
			//std::cout << "imp q, t, x, y, z: " << imp.get_charge() << ", "
			//	<< imp.get_t() << ", " << imp.get_x() << ", " << imp.get_y() 
			//	<< ", " << imp.get_z() << '\n';

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
			ioniz_recomb(imp, bkg, oa_ioniz, oa_recomb, imp_time_step, tidx, 
				xidx, yidx, zidx, ioniz_warnings, recomb_warnings);
		}

	}

	void print_ioniz_recomb_warn(int ioniz_warnings, int recomb_warnings)
	{
		// Number of ionization warnings
		if (ioniz_warnings > 0)
		{
			std::cout << "Warning! The ionization probability was greater than"
				<< " 1.0 for " << ioniz_warnings << " time steps. Consider "
				<< "using a smaller time step.\n";
		}

		// Number of recombination warnings
		if (recomb_warnings > 0)
		{
			std::cout << "Warning! The recombination probability was greater "
				<< "than 1.0 for " << recomb_warnings << " time steps. Consider"
				<< " using a smaller time step.\n";
		}
	}

	void main_loop(const Background::Background& bkg, Statistics& imp_stats,
		const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADAS& oa_recomb)
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
		// time/location to the end. Dynamic scheduling likely the best here 
		// since the loop times can vary widely.
		int thread_imp_count {};
		int thread_imp_num {};
		int ioniz_warnings {};
		int recomb_warnings {};
		#pragma omp parallel for schedule(dynamic) \
			shared(bkg) \
			firstprivate(thread_imp_num, thread_imp_count) \
			reduction(+: imp_stats, ioniz_warnings, recomb_warnings)
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

				// Only do this if enough impurities are being followed, 
				// otherwise the following modulo in the if statement will
				// be a divide by zero error.
				if (thread_imp_num > prog_interval)
				{
					if ((thread_imp_count % (thread_imp_num / prog_interval)) 
						== 0 && thread_imp_count > 0)
					{
						std::cout << "Followed " << thread_imp_count 
							* num_threads << "/" << imp_num << " impurities (" 
							<< static_cast<int>(perc_complete) << "%)\n";
					}
				}
			}

			// Create starting impurity ion
			Impurity primary_imp = create_primary_imp(bkg);
			
			// Each thread starts with one impurity ion, but it may end up 
			// following more than that because that ion may split off another
			// one, which may split off another one, so the thread will
			// continue until there are no more to follow. One could argue that
			// the algorithm here should use OpenMP tasks, which would work 
			// fine, but it is probably unneccesary because all we really need
			// is all the CPUs to stay busy, which is accomplished with the
			// current algorithm with minimal OpenMP overhead. Just some
			// threads may track more impurities than others, which is fine
			// as long as you use dynamic scheduling. 
			std::vector<Impurity> imps {primary_imp};
			while (!imps.empty())
			{
				Impurity imp {pop_back_remove(imps)};

				// Logic for handling additional split impurities not here yet
				follow_impurity(imp, bkg, imp_stats, oa_ioniz, oa_recomb,
					ioniz_warnings, recomb_warnings);
			}

			// Increment thread-specific counter. This just counts primary
			// impurities. 
			thread_imp_count += 1;
		}

		// Let user know how many, if any, warnings occured.
		print_ioniz_recomb_warn(ioniz_warnings, recomb_warnings);
	}

	Statistics follow_impurities(Background::Background& bkg)
	{
		// Initialize particle statistics vectors, all contained within a
		// Statistics object.
		Statistics imp_stats {bkg.get_dim1(), bkg.get_dim2(), 
			bkg.get_dim3(), bkg.get_dim4()};

		// Load OpenADAS data needed for ionization/recombination rates. 
		std::string_view openadas_root {Input::get_opt_str(
			Input::openadas_root)};
		int imp_atom_num {Input::get_opt_int(Input::imp_atom_num)};
		int openadas_year {Input::get_opt_int(Input::openadas_year)};
		OpenADAS::OpenADAS oa_ioniz {openadas_root, openadas_year, 
			imp_atom_num, "scd"};
		OpenADAS::OpenADAS oa_recomb {openadas_root, openadas_year, 
			imp_atom_num, "acd"};

		// Execute main particle following loop.
		main_loop(bkg, imp_stats, oa_ioniz, oa_recomb);
	
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
