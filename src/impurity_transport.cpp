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
#include <cstdlib>

#include "impurity_transport.h"
#include "background.h"
#include "read_input.h"
#include "vectors.h"
#include "impurity.h"
#include "impurity_stats.h"
#include "random.h"
#include "constants.h"
#include "openadas.h"
#include "collisions.h"


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

	double get_birth_y(const Background::Background& bkg, 
		const int imp_ystart_opt_int)
	{
		// Start at specific point
		if (imp_ystart_opt_int == 0)
		{
			return {Input::get_opt_dbl(Input::imp_ystart_val)};
		}

		// Start between a range (here defaults to the full y-width of the 
		// simulation volume).
		else if (imp_ystart_opt_int == 1)
		{
			return Random::get(bkg.get_y_min(), bkg.get_y_max());
		}
		
		return 0.0;
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

	Impurity create_primary_imp(const Background::Background& bkg, 
		const int imp_ystart_opt_int)
	{
		// Load input options as local variables for cleaner code.
		double imp_mass_amu {Input::get_opt_dbl(Input::imp_mass_amu)};
		int imp_atom_num {Input::get_opt_int(Input::imp_atom_num)};

		// Starting t,x,y,z for the impurity
		double t_imp = get_birth_t(bkg);
		double x_imp = get_birth_x(bkg);
		double y_imp = get_birth_y(bkg, imp_ystart_opt_int);
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

	// Delete eventually
	void lorentz_update(Impurity& imp, const Background::Background& bkg,
		const double imp_time_step, const int tidx, const int xidx, 
		const int yidx, const int zidx)
	{
		// Impurity's charge
		double imp_q {imp.get_charge() * -Constants::charge_e};

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
	}

	std::tuple<double, double, double> lorentz_forces(Impurity& imp, 
		const Background::Background& bkg, const int tidx, const int xidx, 
		const int yidx, const int zidx)
	{
		// Impurity's charge
		double imp_q {imp.get_charge() * -Constants::charge_e};

		// Each component of the Lorentz force
		double fx {imp_q * (bkg.get_ex()(tidx, xidx, yidx, zidx) + 
			imp.get_vy() * bkg.get_b()(tidx, xidx, yidx, zidx))};
		double fy {imp_q * (bkg.get_ey()(tidx, xidx, yidx, zidx) + 
			-imp.get_vx() * bkg.get_b()(tidx, xidx, yidx, zidx))};
		double fz {imp_q * bkg.get_ez()(tidx, xidx, yidx, zidx)};

		// Return as tuple
		return std::make_tuple(fx, fy, fz);
	}

	double get_var_time_step_trans(Impurity& imp, 
		const Background::Background& bkg, const int xidx, const int yidx, 
		const int zidx, const double fx, const double fy, const double fz)
	{
		// v0 = velocity before forces are applied
		// v1 = v0 + dv = velocity after forces are applied
		//
		// We use the rule that a particle's step must be less than the
		// minimum width of the cell (w), e.g.,
		// (v0 + dv) * dt < w. We modify this condition to
		// (|v0| + |dv|) * dt < w to avoid imaginary solutions. This in 
		// fact adds an additional layer of conservatism since 
		// |v0 + dv| < |v0| + |dv|, i.e., we are always assuming the
		// "worst case scenario" in which the sign of dv is the same as 
		// v0. This means sometimes we are restricting the time step using
		// too large a velocity, which just means the time step will be
		// smaller than necessary sometimes. This is okay.
		// We add another layer of safety onto this by defining the step 
		// size as a fraction of the cell width, thus,
		//   (|v0| + |dv|) * dt = w / safety_frac    (1)
		//   dv = dt * F / m  (from F = m * dv/dt)  (2)
		// (Note: Using |v0| + |dv| above implies we also should use |F|)
		// Plug (2) into (1) and you get a quadratic equation in dt:
		//   dt^2 + (v0 * m / F)dt - m * w / (safety_frac * F) = 0
		// Quadratic equation says:
		//  dt = (-v0 * m +|- sqrt((v0 * m)^2 
		//		+ 4 * m * w * F / safety_frac))) / (2 * F)
		// If F > 0, we can have two solutions. Common sense says just 
		// pick the one guarunteed to be positive (choose the + in +|-). 
		// If F < 0, then we could get an imaginary solution and be in 
		// trouble. 

		// Fraction of cell width to limit step size to.
		const double safety_frac {5.0};

		// Width of cell in each dimension.
		double xwidth {bkg.get_grid_x()[xidx+1] - bkg.get_grid_x()[xidx]};
		double ywidth {bkg.get_grid_y()[yidx+1] - bkg.get_grid_y()[yidx]};
		double zwidth {bkg.get_grid_z()[zidx+1] - bkg.get_grid_z()[zidx]};
		double min_width {std::min({xwidth, ywidth, zwidth})};
		
		// Velocity and force magnitudes
		double v {std::sqrt(imp.get_vx()*imp.get_vx() 
			+ imp.get_vy()*imp.get_vy() + imp.get_vz()*imp.get_vz())};  
		double f {std::sqrt(fx*fx + fy*fy + fz*fz)};

		// The force can be zero if the particle has recombined to a
		// neutral. This would cause a NaN down below, so we intercept
		// here to calculate the variable time step.
		if (std::abs(f) < Constants::small)
		{
			return min_width / v / safety_frac;
		}

		// The discriminant in the solution for dt
		double disc {v*v * imp.get_mass()*imp.get_mass() 
			+ 4.0 * f * imp.get_mass() * min_width / safety_frac};

		// Check for imaginary solutions
		if (disc < 0)
		{
			std::cerr << "Error! Discriminant in variable time step " 
				<< "calculation is negative and causing an imaginary "
				<< "solution. disc = " << disc << '\n';
		}

		// Two possible solutions. Return the one gauranteed to be 
		// positive (+ disc vs. - disc)?
		double dt {(-v * imp.get_mass() + std::sqrt(disc)) / (2.0 * f)};
		if (std::isnan(dt) || std::isnan(disc) || std::isnan(f) 
			|| std::isnan(v))
		{
			std::cerr << "Error! NaN encountereded in variable time " 
				<< "step.\n" << "v, f, disc, dt = " << v << ", " << f
				<< ", " << disc << ", " << dt << '\n' << "Setting to"
				<< " imp_time_step_min\n";
			std::cerr << "imp_q = " << imp.get_charge() << '\n';
			std::exit(0);
			dt = Input::get_opt_dbl(Input::imp_time_step_min);
		}
		return dt;
	}


	double get_var_time_step(Impurity& imp, 
		const Background::Background& bkg, const int tidx, 
		const int xidx, const int yidx, const int zidx, const double fx, 
		const double fy, const double fz, const bool imp_coll_on)
	{
		// If impurity is at rest, just choose a very low number to kick
		// things off.
		double imp_v {std::sqrt(imp.get_vx()*imp.get_vx() 
			+ imp.get_vy()*imp.get_vy() + imp.get_vz()*imp.get_vz())};
		if (imp_v < Constants::small)
		{
			return 1e-15;
		}
		else
		{
			// Calculate time step for reasonable transport calculation
			double dt_trans {get_var_time_step_trans(imp, bkg, xidx, yidx, 
				zidx, fx, fy, fz)};

			// Calculate time step for a reasonable collisions calculation
			// (if necessary). Don't do this if the impurity is neutral, it
			// won't work. Unfortunately this is calculating the momentum
			// loss factor, tossing it, and then elsehwere we have to 
			// calculate it again. This could use some restructuring to save
			// some time but I'm favoring readibility at this point in time. 
			if (imp_coll_on && imp.get_charge() > 0) 
			{
				double dt_coll {dt_trans};
				Collisions::set_var_time_step_coll(dt_coll, imp, bkg, tidx, 
					xidx, yidx, zidx);

				// Choose the smallest of the two
				//std::cout << "dt_trans, dt_coll = " << dt_trans << ", " 
				//	<< dt_coll << '\n';
				return std::min({dt_trans, dt_coll});
			}

			return dt_trans;
		}
	}

	void step(Impurity& imp, const double fx, const double fy, const double fz, 
		const double imp_time_step)
	{
		// Change in velocity over time step (this is just F = m * dv/dt)
		double dvx {fx * imp_time_step / imp.get_mass()};
		double dvy {fy * imp_time_step / imp.get_mass()};
		double dvz {fz * imp_time_step / imp.get_mass()};

		// Update particle velocity
		imp.set_vx(imp.get_vx() + dvx);
		imp.set_vy(imp.get_vy() + dvy);
		imp.set_vz(imp.get_vz() + dvz);

		// Update particle time and position
		imp.set_t(imp.get_t() + imp_time_step);
		imp.set_x(imp.get_x() + imp.get_vx() * imp_time_step);
		imp.set_y(imp.get_y() + imp.get_vy() * imp_time_step);
		imp.set_z(imp.get_z() + imp.get_vz() * imp_time_step);
	}

	void record_stats(Statistics& imp_stats, const Impurity& imp, 
		const Background::Background& bkg, const int tidx, const int xidx, 
		const int yidx, const int zidx, const double imp_time_step)
	{
		// Add one to counts to this location	
		imp_stats.add_counts(tidx, xidx, yidx, zidx, 1);

		// Add total weight to this location. The total weight is the particle
		// weight times the time step. This is Monte Carlo stuff, best read
		// the literature if you don't recognize it.
		imp_stats.add_weights(tidx, xidx, yidx, zidx, 
			imp.get_weight() * imp_time_step);

		// Add each velocity component to the running sum for this location
		if (imp_stats.get_vel_stats())
		{
			imp_stats.add_vels(tidx, xidx, yidx, zidx, imp.get_vx(), 
				imp.get_vy(), imp.get_vz());
		}

		// Add value of gyroradius to running sum at this location
		imp_stats.add_gyrorad(tidx, xidx, yidx, zidx, imp, bkg);
	}

	bool check_boundary(const Background::Background& bkg, Impurity& imp)
	{
		// If we've run past the time range covered by the background plasma
		// then we're done
		if (imp.get_t() > bkg.get_t_max()) return false; 

		// The x boundaries are treated as absorbing. Could certianly add
		// different options here in the future. There can be Gkeyll
		// artifact near the x bounds, so there is an additional option to 
		// consider anything within imp_xbound_buffer as absorbed.
		const double imp_xbound_buffer {Input::get_opt_dbl(
			Input::imp_xbound_buffer)};
		if (imp.get_x() <= bkg.get_grid_x()[0] + imp_xbound_buffer)
		{
			return false;
		}
		if (imp.get_x() >= bkg.get_grid_x().back() - imp_xbound_buffer) 
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

	void collision(Impurity& imp, const Background::Background& bkg, 
		const double imp_time_step, const int tidx, const int xidx, 
		const int yidx, const int zidx)
	{
		// Local plasma propoerties
		double local_ne = bkg.get_ne()(tidx, xidx, yidx, zidx);
		double local_te = bkg.get_te()(tidx, xidx, yidx, zidx);
		double local_ti = bkg.get_ti()(tidx, xidx, yidx, zidx);

		// Impurity modified within collision step
		Collisions::collision_update(imp, local_te, local_ti, local_ne, 
			imp_time_step);
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
		int& recomb_warnings, const bool imp_coll_on, 
		const bool imp_iz_recomb_on, const int imp_time_step_opt_int)
	{
		// Timestep of impurity transport simulation. It can be a constant
		// value (set here), or set on the fly based on a reasonable criteria.
		// Initialize it with whatever the minimum time step is.
		double imp_time_step {Input::get_opt_dbl(Input::imp_time_step_min)};
		if (imp_time_step_opt_int == 0)
		{
			imp_time_step = Input::get_opt_dbl(Input::imp_time_step);
		}

		// For debugging purposes (only works with one thread)
		//static int imp_id {0};
		//imp_id++;

		bool continue_following {true};
		while (continue_following)
		{
			// Get nearest time index
			int tidx {get_nearest_index(bkg.get_times(), imp.get_t())};

			// Get x,y,z indices
			int xidx {get_nearest_cell_index(bkg.get_grid_x(), imp.get_x())};
			int yidx {get_nearest_cell_index(bkg.get_grid_y(), imp.get_y())};
			int zidx {get_nearest_cell_index(bkg.get_grid_z(), imp.get_z())};

			// Calculate Lorentz force components. First loop these are all
			// zero if particles start at rest.
			auto [fx, fy, fz] = lorentz_forces(imp, bkg, tidx, xidx, yidx, 
				zidx);

			// Calculate variable time step (if necessary)
			if (imp_time_step_opt_int == 1) imp_time_step = 
				get_var_time_step(imp, bkg, tidx, xidx, yidx, zidx, fx, 
				fy, fz, imp_coll_on);

			// Check for a collision
			if (imp_coll_on)
			{
				collision(imp, bkg, imp_time_step, tidx, xidx, yidx, zidx);
			}

			// Debugging
			//std::cout << "id, q, t, x, y, z, dt, fx, fy, fz: " 
			//	<< imp_id << ", " 
			//	<< imp.get_charge() << ", "<< imp.get_t() << ", " 
			//	<< ", " << imp.get_x() << ", " << imp.get_y() 
			//	<< ", " << imp.get_z() << ", " << imp_time_step 
			//	<< ", " << fx << ", " << fy << ", " << fz << '\n';

			// Update statistics. Need to do this after the time step is
			// calculated (if it is), but before the particle moves into 
			// a different cell in step.
			record_stats(imp_stats, imp, bkg, tidx, xidx, yidx, zidx, 
				imp_time_step);

			// Last thing is move particle to a new location
			//step(imp, imp_time_step);
			step(imp, fx, fy, fz, imp_time_step);

			// Check for ionization or recombination
			if (imp_iz_recomb_on)
			{
				ioniz_recomb(imp, bkg, oa_ioniz, oa_recomb, imp_time_step, 
					tidx, xidx, yidx, zidx, ioniz_warnings, recomb_warnings);
			}

			// Check for a terminating or boundary conditions
			continue_following = check_boundary(bkg, imp);
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
		std::string imp_collisions {Input::get_opt_str(Input::imp_collisions)};
		std::string imp_iz_recomb {Input::get_opt_str(Input::imp_iz_recomb)};
		std::string imp_ystart_opt {Input::get_opt_str(Input::imp_ystart_opt)};
		std::string imp_time_step_opt {
			Input::get_opt_str(Input::imp_time_step_opt)};

		// Use internal control variables for string input options. It's better
		// to pass these as booleans or integers to avoid repeatedly checking
		// the value of a string. Should eventually wrap these into a class so
		// they can all be passed via a single reference.
		// 1. Boolean for collisions
		bool imp_coll_on {false};
		if (imp_collisions == "yes") imp_coll_on = true;

		// 2. Boolean for ionization/recombination
		bool imp_iz_recomb_on {false};
		if (imp_iz_recomb == "yes") imp_iz_recomb_on = true;

		// 3. Integer for y start option
		int imp_ystart_opt_int {0};
		if (imp_ystart_opt == "single_value") imp_ystart_opt_int = 0;
		else if (imp_ystart_opt == "range") imp_ystart_opt_int = 1;
		else std::cerr << "Error: Unrecognized value for imp_ystart_opt (" << 
			imp_ystart_opt << "). Valid options are: 'single_value' or " <<
			"'range'.\n";

		// 4. Integer for time step option
		int imp_time_step_opt_int {0};
		if (imp_time_step_opt == "constant") imp_time_step_opt_int = 0;
		else if (imp_time_step_opt == "variable") imp_time_step_opt_int = 1;
		else std::cerr << "Error: Unrecognized value for imp_time_step_opt (" 
			<< imp_time_step_opt << "). Valid options are: 'constant' or " <<
			"'variable'.\n";

		// https://stackoverflow.com/questions/29633531/user-defined-
		// reduction-on-vector-of-varying-size/29660244#29660244
		#pragma omp declare reduction(+: Statistics: \
			omp_out = omp_out + omp_in) \
			initializer(omp_priv(omp_orig))

		// Startup message. Create parallel region to see how many threads it
		// will generate in the following loop.
		std::cout << "Starting particle following...\n";
		#pragma omp parallel
		{
			if (omp_get_thread_num() == 0) std::cout << " Number of threads: " 
				<< omp_get_num_threads() << '\n';
		}

		// Print progress this many times
		constexpr int prog_interval {10};
	
		// Loop through one impurity at a time, tracking it from its birth
		// time/location to the end. Dynamic scheduling likely the best here 
		// since the loop times can vary widely. Declaring bkg, oa_ioniz and
		// oa_recomb as shared is redundant, but I like being explicit. 
		int ioniz_warnings {};
		int recomb_warnings {};
		int imp_count {};
		#pragma omp parallel for schedule(dynamic) \
			shared(bkg, oa_ioniz, oa_recomb) \
			reduction(+: imp_stats, ioniz_warnings, recomb_warnings)
		for (int i = 0; i < imp_num; ++i)
		{

			// Create starting impurity ion
			Impurity primary_imp = create_primary_imp(bkg, imp_ystart_opt_int);
			
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
					ioniz_warnings, recomb_warnings, imp_coll_on, 
					imp_iz_recomb_on, imp_time_step_opt_int);
			}

			// Print out progress at intervals set by prog_interval (i.e.,
			// prog_interval = 10 means print out every 10%, 5 means every 20%,
			// etc.). Would prefer not to have a critical section here, but
			// this is such a quick block of code that it should have very
			// little impact on the overall simulation time considering all
			// the calculation time is spent following impurities. 
			#pragma omp critical
			{
				// Increment shared counter
				++imp_count;
				
				if (imp_num > prog_interval)
				{
					if ((imp_count % (imp_num / prog_interval)) == 0 
						&& imp_count > 0)
					{
						double perc_complete {static_cast<double>(imp_count) 
							/ imp_num * 100};
						std::cout << "Followed " << imp_count << "/" << imp_num 
							<< " impurities (" << static_cast<int>(perc_complete) 
							<< "%)\n";
					}
				}
			}
		}

		// Let user know how many, if any, warnings occured.
		print_ioniz_recomb_warn(ioniz_warnings, recomb_warnings);
	}

	Statistics follow_impurities(Background::Background& bkg)
	{
		// Initialize particle statistics vectors, all contained within a
		// Statistics object. Option to control if the three velocity arrays
		// are allocated (to save memory).
		bool imp_vel_stats {false};
		std::string imp_vel_stats_str 
			{Input::get_opt_str(Input::imp_vel_stats)};
		if (imp_vel_stats_str == "yes") imp_vel_stats = true;
		Statistics imp_stats {bkg.get_dim1(), bkg.get_dim2(), 
			bkg.get_dim3(), bkg.get_dim4(), imp_vel_stats};

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
		// the density by the scale factor.
		std::cout << "Calculating derived quantities...\n";
		int imp_num {Input::get_opt_int(Input::imp_num)};
		double imp_source_scale_fact {
			Input::get_opt_dbl(Input::imp_source_scale_fact)};
		std::cout << "  Density...\n";
		imp_stats.calc_density(bkg, imp_num, imp_source_scale_fact);
		if (imp_stats.get_vel_stats())
		{
			std::cout << "  Velocity...\n";
			imp_stats.calc_vels();
		}
		imp_stats.calc_gyrorad();
		
		return imp_stats;

	}

}
