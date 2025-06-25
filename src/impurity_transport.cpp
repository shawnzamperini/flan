/**
* @file impurity_transport.cpp
*
* @brief Impurity transport routines
*/

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <tuple>
#include <vector>

#include "background.h"
#include "boris.h"
#include "collisions.h"
#include "constants.h"
#include "flan_types.h"
#include "impurity.h"
#include "impurity_stats.h"
#include "impurity_transport.h"
#include "openadas.h"
#include "options.h"
#include "random.h"
#include "variance_reduction.h"
#include "vectors.h"


namespace Impurity
{
	double get_birth_t(const Background::Background& bkg)
	{
		// The random get uses doubles, so need to make sure we are passing
		// double if float is used for BkgFPType
		return Random::get(static_cast<double>(bkg.get_t_min()), 
			static_cast<double>(bkg.get_t_max()));
	}

	double get_birth_x(const Background::Background& bkg,
		const Options::Options& opts)
	{	
		// Uniformily distributed between xmin, xmax.
		double start_x {Random::get(opts.imp_xmin(), opts.imp_xmax())};

		// We need to start the impurity at a grid node, otherwise the 
		// code may autolocate it somewhere further away. This is actually
		// super subtle yet extremely important and can have surprisingly
		// huge and confusing implication if you neglect it. I lost years off 
		// my life figuring this out.
		//int xidx {get_nearest_cell_index(bkg.get_grid_x(), 
		//	static_cast<BkgFPType>(start_x))};
		//return static_cast<double>(bkg.get_grid_x()[xidx]);

		return start_x;
	}

	double get_birth_y(const Background::Background& bkg, 
		const Options::Options& opts)
	{
		// Start at specific point
		double start_y {};
		if (opts.imp_ystart_opt_int() == 0)
		{
			 start_y = opts.imp_ystart_val();
		}

		// Start between a range (here defaults to the full y-width of the 
		// simulation volume).
		else if (opts.imp_ystart_opt_int() == 1)
		{
			start_y = Random::get(static_cast<double>(bkg.get_y_min()), 
				static_cast<double>(bkg.get_y_max()));
		}
		
		// We need to start the impurity at a grid node, otherwise the 
		// code may autolocate it somewhere further away. This is actually
		// super subtle yet extremely important and can have surprisingly
		// huge and confusing implication if you neglect it. I lost years off 
		// my life figuring this out.
		//int yidx {get_nearest_cell_index(bkg.get_grid_y(), 
		//	static_cast<BkgFPType>(start_y))};
		//return static_cast<double>(bkg.get_grid_y()[yidx]);

		return start_y;
	}

	double get_birth_z(const Background::Background& bkg,
		const Options::Options& opts)
	{
		double start_z {};
		if (opts.imp_zstart_opt_int() == 0)
		{
			 start_z = opts.imp_zstart_val();
		}

		// Start between a range (here defaults to the full z-width of the 
		// simulation volume).
		else if (opts.imp_zstart_opt_int() == 1)
		{
			start_z = Random::get(static_cast<double>(bkg.get_z_min()), 
				static_cast<double>(bkg.get_z_max()));
		}

		// We need to start the impurity at a grid node, otherwise the 
		// code may autolocate it somewhere further away. This is actually
		// super subtle yet extremely important and can have surprisingly
		// huge and confusing implication if you neglect it. I lost years off 
		// my life figuring this out.
		//int zidx {get_nearest_cell_index(bkg.get_grid_z(), 
		//	static_cast<BkgFPType>(start_z))};
		//return static_cast<double>(bkg.get_grid_z()[zidx]);

		return start_z;
	}

	int get_birth_charge(const Options::Options& opts)
	{
		return opts.imp_init_charge();
	}

	Impurity create_primary_imp(const Background::Background& bkg, 
		const Options::Options& opts)
	{
		// Starting t,x,y,z for the impurity
		double t_imp = get_birth_t(bkg);
		double x_imp = get_birth_x(bkg, opts);
		double y_imp = get_birth_y(bkg, opts);
		double z_imp = get_birth_z(bkg, opts);
		auto [X_imp, Y_imp, Z_imp] = opts.mapc2p()(x_imp, y_imp, z_imp);

		// Assume starting at rest, but if this ever changes do it here
		double vX_imp {0.0};
		double vY_imp {0.0};
		double vZ_imp {0.0};

		// Impurity starting charge
		int charge_imp = get_birth_charge(opts);
	
		// Return a temporary Impurity object. Since these are primary Impurity
		// objects they by definition start with weight = 1.0.
		double imp_weight {1.0};
		double imp_mass {opts.imp_mass_amu() * Constants::amu_to_kg};
		return {t_imp, x_imp, y_imp, z_imp, X_imp, Y_imp, Z_imp, vX_imp, 
			vY_imp, vZ_imp, imp_weight, charge_imp, imp_mass, 
			opts.imp_atom_num()};
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
		// In this example, we want 1 returned, so we return 2 - 1 = 1. If
		// value is outside the range, return the index of the respective end
		// of the vector.
		int index = std::distance(grid_edges.begin(), lower);

		// If less than everything, return the first cell.
		if (lower == grid_edges.begin()) return 0;

		// If larger than everything, return the last cell. Note end() is the
		// iterator that points past the last element of a vector, so to return
		// the cell we need to subtract by 2. 
		else if (lower == grid_edges.end()) 
		{
			return index - 2;
		}

		//else return lower - grid_edges.begin() - 1;
		else return index - 1;
	}

	std::tuple<double, double, double> lorentz_forces(Impurity& imp, 
		const Background::Background& bkg, const int tidx, const int xidx, 
		const int yidx, const int zidx)
	{
		// Impurity's charge
		double imp_q {imp.get_charge() * -Constants::charge_e};

		// Electric and magnetic field components
		double eX {bkg.get_eX()(tidx, xidx, yidx, zidx)};
		double eY {bkg.get_eY()(tidx, xidx, yidx, zidx)};
		double eZ {bkg.get_eZ()(tidx, xidx, yidx, zidx)};
		double bX {bkg.get_bX()(tidx, xidx, yidx, zidx)};
		double bY {bkg.get_bY()(tidx, xidx, yidx, zidx)};
		double bZ {bkg.get_bZ()(tidx, xidx, yidx, zidx)};

		// Impurity velocity components
		double vX {imp.get_vX()};
		double vY {imp.get_vY()};
		double vZ {imp.get_vZ()};

		//std::cout << "  Ex, Ey, Ez = " << eX << ", " << eY << ", " << eZ 
		//	<< "\n";
		//std::cout << "  Bx, By, Bz, B = " << bX << ", " << bY << ", " << bZ 
		//	<< ", " << std::sqrt(bX*bX + bY*bY + bZ*bZ) << "\n";
		//std::cout << "  vX, vY, vZ = " << vX << ", " << vY << ", " << vZ 
		//	<< "\n";

		// Each component of the Lorentz force in physical space
		double fX {imp_q * (eX + vY * bZ - vZ * bY)};
		double fY {imp_q * (eY + vZ * bX - vX * bZ)};
		double fZ {imp_q * (eZ + vX * bY - vY * bX)};

		// Return as tuple
		return std::make_tuple(fX, fY, fZ);
	}

	double get_var_time_step_trans(Impurity& imp, 
		const Background::Background& bkg, const int xidx, const int yidx, 
		const int zidx, const double fx, const double fy, const double fz,
		const Options::Options& opts)
	{
		return opts.imp_time_step();
		/*
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
		
		// Don't forget to sqrt it.
		min_width = std::sqrt(min_width);

		// Velocity and force magnitudes
		double v {std::sqrt(imp.get_vX()*imp.get_vX() 
			+ imp.get_vY()*imp.get_vY() + imp.get_vZ()*imp.get_vZ())};  
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
			dt = opts.imp_time_step_min();
		}
		return dt;
		*/
	}

	double get_var_time_step(Impurity& imp, 
		const Background::Background& bkg, const int tidx, 
		const int xidx, const int yidx, const int zidx, const double fx, 
		const double fy, const double fz, const Options::Options& opts)
	{
		// A temporary thing
		//std::cerr << "Error! Variable time step is not working yet with new"
		//	<< " geometry. Defaulting to input value.\n";
		//return opts.imp_time_step();

		// If impurity is at rest, just choose a very low number to kick
		// things off.
		double imp_v {std::sqrt(imp.get_vX()*imp.get_vX() 
			+ imp.get_vY()*imp.get_vY() + imp.get_vZ()*imp.get_vZ())};
		if (imp_v < Constants::small)
		{
			return 1e-15;
		}
		else
		{
			// Calculate time step for reasonable transport calculation
			//double dt_trans {get_var_time_step_trans(imp, bkg, xidx, yidx, 
			//	zidx, fx, fy, fz, opts)};

			// Unsure how to calculate this, so just assign to whatever is 
			// input.
			double dt_trans {opts.imp_time_step()};

			// Calculate time step for a reasonable collisions calculation
			// (if necessary). Don't do this if the impurity is neutral, it
			// won't work. Unfortunately this is calculating the momentum
			// loss factor, tossing it, and then elsewhere we have to 
			// calculate it again. This could use some restructuring to save
			// some time but I'm favoring readibility at this point in time. 
			if (opts.imp_collisions_int() > 0 && imp.get_charge() > 0) 
			{
				double dt_coll {dt_trans};
				Collisions::set_var_time_step_coll(dt_coll, imp, bkg, tidx, 
					xidx, yidx, zidx, opts);

				// Choose the smallest of the two
				//std::cout << "dt_trans, dt_coll = " << dt_trans << ", " 
				//	<< dt_coll << '\n';
				return std::min({dt_trans, dt_coll});
			}

			return dt_trans;
		}
	}

	bool step(Impurity& imp, const double fX, const double fY, const double fZ, 
		const double imp_time_step, const Background::Background& bkg,
		const Options::Options& opts, const int tidx, int& xidx, int& yidx, 
		int& zidx)
	{
		// If we've run past the time range covered by the background plasma
		// then we're done
		imp.set_t(imp.get_t() + imp_time_step);
		if (imp.get_t() > bkg.get_t_max())
		{
			return false;
		}

		// Change in velocity over time step (this is just F = m * dv/dt)
		//double dvX {fX * imp_time_step / imp.get_mass()};
		//double dvY {fY * imp_time_step / imp.get_mass()};
		//double dvZ {fZ * imp_time_step / imp.get_mass()};

		// Update particle velocity
		//imp.set_vX(imp.get_vX() + dvX);
		//imp.set_vY(imp.get_vY() + dvY);
		//imp.set_vZ(imp.get_vZ() + dvZ);

		// Use Boris algorithm to update particle velocity
		Boris::update_velocity(imp, bkg, imp_time_step, tidx, xidx, yidx, 
			zidx);

		// Update particle position in physical space.
		imp.set_X(imp.get_X() + imp.get_vX() * imp_time_step);
		imp.set_Y(imp.get_Y() + imp.get_vY() * imp_time_step);
		imp.set_Z(imp.get_Z() + imp.get_vZ() * imp_time_step);

        // Using the metric coefficients, update the particle position in 
        // computational space. The equation is derived from Dhaeseleer 2.4.2.
        // The equations are:
        //   dx = g00*dX + g10*dY + g20*dZ
        //   dy = g01*dX + g11*dY + g21*dZ
        //   dz = g02*dX + g12*dY + g22*dZ
        // Where I've already accounted for the fact that gij is symmetric and
        // swapped variables according (e.g., gij01 = gij10). 
        double dX {imp.get_X() - imp.get_prevX()};
        double dY {imp.get_Y() - imp.get_prevY()};
        double dZ {imp.get_Z() - imp.get_prevZ()};
        double dx {bkg.get_gij_00()(xidx, yidx, zidx) * dX
            + bkg.get_gij_01()(xidx, yidx, zidx) * dY 
            + bkg.get_gij_02()(xidx, yidx, zidx) * dZ};
        double dy {bkg.get_gij_01()(xidx, yidx, zidx) * dX
            + bkg.get_gij_11()(xidx, yidx, zidx) * dY 
            + bkg.get_gij_22()(xidx, yidx, zidx) * dZ};
        double dz {bkg.get_gij_02()(xidx, yidx, zidx) * dX
            + bkg.get_gij_12()(xidx, yidx, zidx) * dY 
            + bkg.get_gij_22()(xidx, yidx, zidx) * dZ};
        imp.set_x(imp.get_x() + dx);
        imp.set_y(imp.get_y() + dy);
        imp.set_z(imp.get_z() + dz);

        // Bound checking (move to separate function). 
		// Absorbing at maximum x
        //if (imp.get_x() < bkg.get_grid_x()[0] || 
        //    imp.get_x() > bkg.get_grid_x().back()) return false;
		if (imp.get_x() > bkg.get_grid_x().back()) return false;

		// At minimum x move the particle to a random y,z cell. This is a
		// rough approximation to entering the core and leaving it 
		// somewhere else.
		else if (imp.get_x() < bkg.get_grid_x()[0])
		{
			imp.set_y(Random::get(static_cast<double>(bkg.get_y_min()), 
				static_cast<double>(bkg.get_y_max())));
			imp.set_z(Random::get(static_cast<double>(bkg.get_z_min()), 
				static_cast<double>(bkg.get_z_max())));
		}

        // Periodic y boundary
        if (imp.get_y() < bkg.get_grid_y()[0])
        {
            imp.set_y(bkg.get_grid_y().back() + (imp.get_y() 
				- bkg.get_grid_y()[0]));
        }
        else if (imp.get_y() > bkg.get_grid_y().back())
        {
            imp.set_y(bkg.get_grid_y()[0] + (imp.get_y() 
				- bkg.get_grid_y().back()));
        }

		// Absorbing z boundary in SOL, periodic in core
		if (imp.get_x() > opts.lcfs_x())
		{
			if (imp.get_z() < bkg.get_grid_z()[0] || 
				imp.get_z() > bkg.get_grid_z().back()) return false;
		}
		else
		{
			if (imp.get_z() < bkg.get_grid_z()[0])
			{
				imp.set_z(bkg.get_grid_z().back() + (imp.get_z() 
					- bkg.get_grid_z()[0]));
			}
			else if (imp.get_z() > bkg.get_grid_z().back())
			{
				imp.set_z(bkg.get_grid_z()[0] + (imp.get_z() 
					- bkg.get_grid_z().back()));
			}
		}

		return true;
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
			static_cast<BkgFPType>(imp.get_weight() * imp_time_step));

		// Add each velocity component to the running sum for this location
		if (imp_stats.get_vel_stats())
		{
			imp_stats.add_vels(tidx, xidx, yidx, zidx, 
				static_cast<BkgFPType>(imp.get_vX()), 
				static_cast<BkgFPType>(imp.get_vY()), 
				static_cast<BkgFPType>(imp.get_vZ()));
		}

		// Add value of gyroradius to running sum at this location
		//imp_stats.add_gyrorad(tidx, xidx, yidx, zidx, imp, bkg);
	}

	void find_containing_cell(Impurity& imp, 
		const Background::Background& bkg, 
		int& xidx, int& yidx, int& zidx, const bool debug)
	{
        xidx = get_nearest_cell_index(bkg.get_grid_x(),     
            static_cast<BkgFPType>(imp.get_x()));
        yidx = get_nearest_cell_index(bkg.get_grid_y(), 
            static_cast<BkgFPType>(imp.get_y()));
        zidx = get_nearest_cell_index(bkg.get_grid_z(), 
            static_cast<BkgFPType>(imp.get_z()));
	}

	void collision(Impurity& imp, const Background::Background& bkg, 
		const double imp_time_step, const int tidx, const int xidx, 
		const int yidx, const int zidx, const Options::Options& opts)
	{
		// Local plasma propoerties
		double local_ne = bkg.get_ne()(tidx, xidx, yidx, zidx);
		double local_te = bkg.get_te()(tidx, xidx, yidx, zidx);
		double local_ti = bkg.get_ti()(tidx, xidx, yidx, zidx);

		// Impurity modified within collision step
		Collisions::collision_update(imp, local_te, local_ti, local_ne, 
			imp_time_step, opts);
	}

	void follow_impurity(Impurity& imp, const Background::Background& bkg, 
		Statistics& imp_stats, const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADAS& oa_recomb, int& ioniz_warnings, 
		int& recomb_warnings, std::vector<Impurity>& imps,
		const std::vector<int> var_red_counts, 
		const bool var_red_on, 
		const Options::Options& opts)
	{
		// Timestep of impurity transport simulation. It can be a constant
		// value (set here), or set on the fly based on a reasonable criteria.
		// Initialize it with whatever the minimum time step is.
		double imp_time_step {opts.imp_time_step_min()};
		if (opts.imp_time_step_opt_int() == 0)
		{
			imp_time_step = opts.imp_time_step();
		}

		// Get nearest time index
		int tidx {get_nearest_index(bkg.get_times(), 
			static_cast<BkgFPType>(imp.get_t()))};

		// Get x,y,z indices in computational space
		int xidx {get_nearest_cell_index(bkg.get_grid_x(), 
			static_cast<BkgFPType>(imp.get_x()))};
		int yidx {get_nearest_cell_index(bkg.get_grid_y(), 
			static_cast<BkgFPType>(imp.get_y()))};
		int zidx {get_nearest_cell_index(bkg.get_grid_z(), 
			static_cast<BkgFPType>(imp.get_z()))};
		find_containing_cell(imp, bkg, xidx, yidx, zidx);

		// Boris algorithm: The velocity stored in the Impurity objects
		// will actually be the velocity at a half timestep earlier, i.e.,
		// at t - dt/2. So we still need to push the particle velocity
		// back by half a time step.
		Boris::update_velocity(imp, bkg, -imp_time_step / 2.0, tidx, xidx, 
			yidx, zidx);

		// Record starting position in statistics arrays
		record_stats(imp_stats, imp, bkg, tidx, xidx, yidx, zidx, 
			imp_time_step);

		// For debugging purposes (only works with one thread)
		static int imp_id {0};
		imp_id++;
		bool debug {false};

		// Begin main particle following loop
		bool continue_following {true};
		while (continue_following)
		{
			// Get nearest time index
			int tidx {get_nearest_index(bkg.get_times(), 
				static_cast<BkgFPType>(imp.get_t()))};

			// Check for ionization or recombination
			if (opts.imp_iz_recomb_int() > 0)
			{
				OpenADAS::ioniz_recomb(imp, bkg, oa_ioniz, oa_recomb, 
					imp_time_step, tidx, xidx, yidx, zidx, ioniz_warnings, 
					recomb_warnings);
			}

			// Calculate Lorentz force components. First loop these are all
			// zero if particles start at rest.
			auto [fX, fY, fZ] = lorentz_forces(imp, bkg, tidx, xidx, yidx, 
				zidx);

			// Variance reduction scheme. This may change the particle's 
			// weight and create a secondary Impurity to be followed. It
			// currently is based on ionization/recombination, hence that also
			// needs to be on to work.
			if (var_red_on && opts.imp_iz_recomb_int() > 0)
			{
				VarianceReduction::split_particle_main(imp, tidx, xidx, yidx, 
					zidx, imp_stats, imps, opts.var_red_min_weight(), 
					var_red_counts, bkg, oa_ioniz, oa_recomb, 
					imp_time_step);
			}

			// Calculate variable time step (if necessary). NOT TRUSTWORTHY
			// YET AFTER UPGRADING PARTICLE FOLLOWING.
			if (opts.imp_time_step_opt_int() == 1) imp_time_step = 
				get_var_time_step(imp, bkg, tidx, xidx, yidx, zidx, fX, 
				fY, fZ, opts);

			// Check for a collision
			if (opts.imp_collisions_int() > 0)
			{
				collision(imp, bkg, imp_time_step, tidx, xidx, yidx, zidx, 
					opts);
			}

			// Debugging
			if (debug)
			{
				double R {std::sqrt(imp.get_X()*imp.get_X() 
					+ imp.get_Y()*imp.get_Y() + imp.get_Z()*imp.get_Z())};
				std::cout << "Before step\n";
				std::cout << "  id, tidx, q, t, x, y, z, dt, fX, fY, fZ, " << 
					"vX, vY, vZ, X, Y, Z, R: \n" << " " << imp_id << ", " 
					<< tidx << ", " << imp.get_charge() << ", "<< imp.get_t() 
					<< ", " << imp.get_x() << ", " << imp.get_y() << ", " 
					<< imp.get_z() << ", " << imp_time_step << '\n' << "  " 
					<< fX << ", " << fY << ", " << fZ << ", " << imp.get_vX() 
					<< ", " << imp.get_vY() << ", " << imp.get_vZ() << ", "
					<< imp.get_X() << ", " << imp.get_Y() << ", " 
					<< imp.get_Z() << ", " << R << '\n' << tidx << " " 
					<< xidx << " " << yidx << " " << zidx << '\n';
			}

			// Perform particle step, update time and tidx. The final location 
			// of the particle could end up out of the grid. That's what the 
			// following code does before recording stats.
			continue_following = step(imp, fX, fY, fZ, imp_time_step, bkg, 
				opts, tidx, xidx, yidx, zidx);

			if (debug)
			{
				double R {std::sqrt(imp.get_X()*imp.get_X() 
					+ imp.get_Y()*imp.get_Y() + imp.get_Z()*imp.get_Z())};
				std::cout << "After step\n";
				std::cout << "  id, tidx, q, t, x, y, z, dt, fX, fY, fZ, " << 
					"vX, vY, vZ, X, Y, Z, R: \n" << " " << imp_id << ", " 
					<< tidx << ", " << imp.get_charge() << ", "<< imp.get_t() 
					<< ", " << imp.get_x() << ", " << imp.get_y() << ", " 
					<< imp.get_z() << ", " << imp_time_step << '\n' << "  " 
					<< fX << ", " << fY << ", " << fZ << ", " << imp.get_vX() 
					<< ", " << imp.get_vY() << ", " << imp.get_vZ() << ", "
					<< imp.get_X() << ", " << imp.get_Y() << ", " 
					<< imp.get_Z() << ", " << R << '\n' << tidx << " " 
					<< xidx << " " << yidx << " " << zidx << '\n';
			}

			// Update xidx, yidx, zidx with the new position
			find_containing_cell(imp, bkg, xidx, yidx, zidx);
			
			// Record stats in the cell the particle ended in. Note: This is 
			// only okay for constant time steps. If it is variable, then cells 
			// that use small time steps could be unfairly influenced by 
			// neighboring cells with larger time steps. We would need to 
			// figure out how much of the particle's step was in each cell and 
			// divide the weight up accordingly. Forgoing this for now, but 
			// it's a future consideration.
			if (continue_following) record_stats(imp_stats, imp, bkg, tidx, 
				xidx, yidx, zidx, imp_time_step);
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
		const OpenADAS::OpenADAS& oa_recomb, 
		const Options::Options& opts)
	{
		// Variance reduction requires ionization/recombination be on
		//if (var_red_on && !imp_iz_recomb_on)
		if (opts.var_red_int() > 0 && opts.imp_iz_recomb_int() == 0)
		{
			std::cout << "Warning! Variance reduction requires "
				<< "imp_ioniz_recomb be set to 'yes'. Variance reduction will"
				<< " not occur for this simulation.\n";
		}

		// Vector of counts at each frame, below which is considered a 
		// low-count region (as defined in variance_reduction.get_counts).
		std::vector<int> var_red_counts (
			imp_stats.get_counts().get_dim1(), 0);

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
		//constexpr int prog_interval {10};
		const int prog_interval {opts.print_interval()};
	
		// Loop through one impurity at a time, tracking it from its birth
		// time/location to the end. Dynamic scheduling likely the best here 
		// since the loop times can vary widely. Declaring bkg, oa_ioniz and
		// oa_recomb as shared is redundant, but I like being explicit. 
		int ioniz_warnings {};
		int recomb_warnings {};
		int prim_imp_count {};
		int tot_imp_count {};
		int priv_count {};
		bool var_red_control {false};
		#pragma omp parallel for schedule(dynamic) \
			firstprivate(priv_count, var_red_counts) \
			firstprivate(var_red_control) \
			shared(bkg, oa_ioniz, oa_recomb) \
			reduction(+: imp_stats, ioniz_warnings, recomb_warnings) \
			//reduction(+: tot_imp_count)
		for (int i = 0; i < opts.imp_num(); ++i)
		{

			// Periodically determine what number of counts qualifies for 
			// splitting a particle, if necessary. 
			if (opts.var_red_int() > 0 && priv_count > 0)
			{
				// This variable is just an integer saying "every X particles
				// followed, update the variance reduction counts". It needs to
				// be greater than zero.
				int var_red_med_mod {static_cast<int>(
					static_cast<double>(opts.imp_num()) / omp_get_num_threads() 
					* opts.var_red_freq())};
				if (var_red_med_mod == 0) var_red_med_mod = 1;

				// Update variance reduction counts according to the variance
				// reduction frequency (var_red_freq), which is a number 
				// between 0 and 1. So for example, a value of 0.10 would 
				// execute this if statement 10 times during the simulation.
				// The modifier determines what fraction of the median is
				// considered low-count. Lower values = more splitting.
				if (priv_count % var_red_med_mod == 0 && priv_count > 0)
				{
					var_red_counts = 
						VarianceReduction::get_counts(imp_stats, 
						var_red_med_mod);

					// We use this separate boolean (initially false) to
					// control if variance reduction occurs so that the
					// simulation can run through for a bit and learn
					// what a low-count region qualifies as. 
					var_red_control = true;

					//for (auto c : var_red_counts)
					//{
					//	std::cout << "var_red_counts = " << c << '\n';
					//}
				}
			}	

			// Create starting impurity ion
			Impurity primary_imp = create_primary_imp(bkg, opts);
			
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
				// Grab an Impurity to follow. First one will be the primary
				// Impurity, any added during the course of follow_impurity
				// are secondary Impuritys split off from the primary.
				Impurity imp {pop_back_remove(imps)};

				// Begin following Impurity
				follow_impurity(imp, bkg, imp_stats, oa_ioniz, oa_recomb,
					ioniz_warnings, recomb_warnings, imps,
					var_red_counts, var_red_control, 
					opts);
				//std::cout << "imps.size() = " << imps.size() << '\n';
				#pragma omp critical
				{
					++tot_imp_count;
				}
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
				++prim_imp_count;
				
				if (opts.imp_num() > prog_interval)
				{
					if ((prim_imp_count % (opts.imp_num() / prog_interval)) 
						== 0 && prim_imp_count > 0)
					{
						double perc_complete {static_cast<double>(
							prim_imp_count) / opts.imp_num() * 100};
						std::cout << "Followed " << prim_imp_count << "/" 
							<< opts.imp_num() << " primary impurities (" 
							<< static_cast<int>(perc_complete) << "%)\n";

						// Give user an update on secondary impurities, since
						// most the time can actually be spent following these.
						std::cout << "  Secondary impurities followed: " 
							<< tot_imp_count - prim_imp_count << '\n';
					}
				}
			}

			// Increment counter private to each thread
			priv_count++;
		}

		// Let user know how many, if any, warnings occured.
		print_ioniz_recomb_warn(ioniz_warnings, recomb_warnings);
	}

	Statistics follow_impurities(Background::Background& bkg, 
		const Options::Options& opts)
	{
		// Initialize particle statistics vectors, all contained within a
		// Statistics object. Option to control if the three velocity arrays
		// are allocated (to save memory).
		Statistics imp_stats {bkg.get_dim1(), bkg.get_dim2(), 
			bkg.get_dim3(), bkg.get_dim4(), 
			static_cast<bool>(opts.imp_vel_stats_int())};

		// Load OpenADAS data needed for ionization/recombination rates. 
		OpenADAS::OpenADAS oa_ioniz {opts.openadas_root(), 
			opts.openadas_year(), opts.imp_atom_num(), "scd"};
		OpenADAS::OpenADAS oa_recomb {opts.openadas_root(), 
			opts.openadas_year(), opts.imp_atom_num(), "acd"};

		// Execute main particle following loop.
		main_loop(bkg, imp_stats, oa_ioniz, oa_recomb, opts);
	
		// Convert the statistics into meaningful quantities. We are scaling
		// the density by the scale factor.
		std::cout << "Calculating derived quantities...\n";
		std::cout << "  Density...\n";
		imp_stats.calc_density(bkg, opts.imp_num(), 
			opts.imp_source_scale_fact());
		if (imp_stats.get_vel_stats())
		{
			std::cout << "  Velocity...\n";
			imp_stats.calc_vels();
		}
		imp_stats.calc_gyrorad();
		
		return imp_stats;

	}

}
