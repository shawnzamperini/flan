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
	// Function to return a random starting time for an impurity ion. This is
	// obviously very basic, to have as its own function, but we do this since
	// this will probably get more complicated in the future.
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

	

	// Create a primary Impurity. Primary impurities are those that start
	// according to the initial condition options specified in the input file.
	// Later, these can create secondary Impurity objects as a result of
	// Monte Carlo splitting.
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

	// Helper function to get the last element in a vector while also removing
	// it.
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
		
		auto lower = std::lower_bound(grid_edges.begin(), grid_edges.end(), 
			value);

		// Comment here explaining this logic again
		if (lower == grid_edges.begin()) return 0;
		else return lower - grid_edges.begin() - 1;
	}

	void do_lorentz_step(Impurity& imp, const Background::Background& bkg,
		const double imp_time_step, const int tidx, const int xidx, 
		const int yidx, const int zidx)
	{
		// Impurity's charge
		double imp_q {imp.get_charge() * Constants::charge_e};

		// Each component of the Lorentz force
		double fx {imp_q * (bkg.get_ex()(tidx, xidx, yidx, zidx) + 
			imp.get_vy() * bkg.get_b()(tidx, xidx, yidx, zidx))};
		double fy {imp_q * (bkg.get_ey()(tidx, xidx, yidx, zidx) + 
			imp.get_vx() * bkg.get_b()(tidx, xidx, yidx, zidx))};
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

		// Add particle's weight to this location
		imp_stats.add_weights(tidx, xidx, yidx, zidx, imp.get_weight());
	}

	bool check_boundary(const Background::Background& bkg, Impurity& imp)
	{
		// The x boundaries are treated as absorbing. Could certianly add
		// different options here in the future.
		if (imp.get_x() <= bkg.get_grid_x()[0]) return false;
		if (imp.get_x() >= bkg.get_grid_x().back()) return false;

		// Absorbing z boundaries
		if (imp.get_z() <= bkg.get_grid_z()[0]) return false;
		if (imp.get_z() >= bkg.get_grid_z().back()) return false;

		// Perioidic y boundaries. This means if an impurity crosses one of
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
			// Get nearest time index
			int tidx {get_nearest_index(bkg.get_times(), imp.get_t())};

			// Get x,y,z indices
			int xidx {get_nearest_cell_index(bkg.get_grid_x(), imp.get_x())};
			int yidx {get_nearest_cell_index(bkg.get_grid_y(), imp.get_y())};
			int zidx {get_nearest_cell_index(bkg.get_grid_z(), imp.get_z())};

			// Perform a step according to the Lorentz force
			do_lorentz_step(imp, bkg, imp_time_step, tidx, xidx, yidx, zidx);

			// Update paeticle time
			imp.set_t(imp.get_t() + imp_time_step);

			// Update statistics
			record_stats(imp_stats, imp, tidx, xidx, yidx, zidx);

			// Check for a terminating or boundary conditions
			continue_following = check_boundary(bkg, imp);

			// Check for a collision
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
		#pragma omp parallel for schedule(dynamic) \
			shared(bkg) \
			reduction(+: imp_stats)
		for (int i = 0; i < imp_num; ++i)
		{
			// OpenMP overhead
			[[maybe_unused]] int thread_id {omp_get_thread_num()};
			[[maybe_unused]] int num_threads {omp_get_num_threads()};

			if (thread_id == 0)
			{
				std::cout << "Number of threads: " << num_threads << '\n';
			}

			// Create starting impurity ion
			Impurity primary_imp = create_primary_imp(bkg);
			
			std::vector<Impurity> imps {primary_imp};
			while (!imps.empty())
			{
				Impurity imp {pop_back_remove(imps)};

				follow_impurity(imp, bkg, imp_stats);
			}
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
	
		// Aggregate the statistics into meaningful quantities.
		
		return imp_stats;

	}

}
