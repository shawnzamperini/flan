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
#include "collisions.h"
#include "constants.h"
#include "impurity.h"
#include "impurity_stats.h"
#include "impurity_transport.h"
#include "kdtree.h"
#include "moeller_trumbore.h"
#include "openadas.h"
#include "options.h"
#include "random.h"
#include "variance_reduction.h"
#include "vectors.h"


namespace Impurity
{
	double get_birth_t(const Background::Background& bkg)
	{
		return Random::get(bkg.get_t_min(), bkg.get_t_max());
	}

	double get_birth_x(const Background::Background& bkg,
		const Options::Options& opts)
	{	
		// Uniformily distributed between xmin, xmax.
		double start_x {Random::get(opts.imp_xmin(), opts.imp_xmax())};

		// We need to start the impurity at the cell center, otherwise the 
		// code my autolocate it somewhere further away.
		int xidx {get_nearest_cell_index(bkg.get_grid_x(), start_x)};
		return bkg.get_x()[xidx];

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
			start_y = Random::get(bkg.get_y_min(), bkg.get_y_max());
		}
		
		// We need to start the impurity at the cell center, otherwise the 
		// code my autolocate it somewhere further away.
		int yidx {get_nearest_cell_index(bkg.get_grid_y(), start_y)};
		return bkg.get_y()[yidx];
	}

	double get_birth_z(const Background::Background& bkg,
		const Options::Options& opts)
	{
		double start_z {opts.imp_zstart_val()};
		int zidx {get_nearest_cell_index(bkg.get_grid_z(), start_z)};
		return bkg.get_z()[zidx];
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
		// In this example, we want 1 returned, so we return 2 - 1 = 1. 
		if (lower == grid_edges.begin()) return 0;
		else return lower - grid_edges.begin() - 1;
	}

	void lorentz_update(Impurity& imp, const Background::Background& bkg,
		const double imp_time_step, const int tidx, const int xidx, 
		const int yidx, const int zidx)
	{
		// Each component of the Lorentz force in physical space
		auto [fX, fY, fZ] = lorentz_forces(imp, bkg, tidx, xidx, yidx, zidx);

		// Change in velocity over time step
		double dvX {fX * imp_time_step / imp.get_mass()};
		double dvY {fY * imp_time_step / imp.get_mass()};
		double dvZ {fZ * imp_time_step / imp.get_mass()};

		// Update particle velocity in physical space
		imp.set_vX(imp.get_vX() + dvX);
		imp.set_vY(imp.get_vY() + dvY);
		imp.set_vZ(imp.get_vZ() + dvZ);
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
		//double fX {imp_q * (bkg.get_eX()(tidx, xidx, yidx, zidx) + 
		//	imp.get_vY() * bkg.get_b()(tidx, xidx, yidx, zidx))
		//	-imp.get_vZ() * bkg.get_b()(tidx, xidx, yidx, zidx)};
		//double fY {imp_q * (bkg.get_eY()(tidx, xidx, yidx, zidx) + 
		//	-imp.get_vX() * bkg.get_b()(tidx, xidx, yidx, zidx))};
		//double fZ {imp_q * bkg.get_eZ()(tidx, xidx, yidx, zidx)};

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
		std::cerr << "Error! Variable time step is not working yet with new"
			<< " geometry. Defaulting to input value.\n";
		return opts.imp_time_step();

		/*
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
				zidx, fx, fy, fz, opts)};

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
		*/
	}

	void step(Impurity& imp, const double fX, const double fY, const double fZ, 
		const double imp_time_step, const Background::Background& bkg,
		std::unique_ptr<KDTree::KDTree_t>& kdtree)
	{
		// Change in velocity over time step (this is just F = m * dv/dt)
		double dvX {fX * imp_time_step / imp.get_mass()};
		double dvY {fY * imp_time_step / imp.get_mass()};
		double dvZ {fZ * imp_time_step / imp.get_mass()};

		// Update particle velocity
		imp.set_vX(imp.get_vX() + dvX);
		imp.set_vY(imp.get_vY() + dvY);
		imp.set_vZ(imp.get_vZ() + dvZ);

		// Update particle time and position in physical space
		imp.set_t(imp.get_t() + imp_time_step);
		imp.set_X(imp.get_X() + imp.get_vX() * imp_time_step);
		imp.set_Y(imp.get_Y() + imp.get_vY() * imp_time_step);
		imp.set_Z(imp.get_Z() + imp.get_vZ() * imp_time_step);

		// Find what index in the physical coordinates is closest
		//std::cout << "Finding nearest neighbors...\n";
		std::size_t nearest_idx {KDTree::nearest_neighbor(kdtree, imp.get_X(),
			imp.get_Y(), imp.get_Z())};

		// Get the computational coordinate indices and update particle
		// position with these values. Can use any of X,Y,Z, just need it for
		// the get_i, etc. function which is the same among each.
		int xidx {bkg.get_X().get_i(nearest_idx)};
		int yidx {bkg.get_X().get_j(nearest_idx)};
		int zidx {bkg.get_X().get_k(nearest_idx)};
		imp.set_x(bkg.get_x()[xidx]);
		imp.set_y(bkg.get_y()[yidx]);
		imp.set_z(bkg.get_z()[zidx]);
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
			imp_stats.add_vels(tidx, xidx, yidx, zidx, imp.get_vX(), 
				imp.get_vY(), imp.get_vZ());
		}

		// Add value of gyroradius to running sum at this location
		//imp_stats.add_gyrorad(tidx, xidx, yidx, zidx, imp, bkg);
	}

	// Get bounding quadrilateral in the x direction in physical coordinates.
	// Must go around the edge of the grid, no cutting across the middle.
	std::array <double, 12> get_x_bound_vertices(
		const Background::Background& bkg, const Options::Options& opts,
		const int xidx, const int yidx, const int zidx)
	{
		// Vertex #1 @ (xidx, yidx, zidx)
		auto [v1x, v1y, v1z] = opts.mapc2p()(bkg.get_grid_x()[xidx], 
			bkg.get_grid_y()[yidx], bkg.get_grid_z()[zidx]);

		// Vertex #2 @ (xidx, yidx+1, zidx)
		auto [v2x, v2y, v2z] = opts.mapc2p()(bkg.get_grid_x()[xidx], 
			bkg.get_grid_y()[yidx+1], bkg.get_grid_z()[zidx]);

		// Vertex #3 @ (xidx, yidx+1, zidx+1)
		auto [v3x, v3y, v3z] = opts.mapc2p()(bkg.get_grid_x()[xidx], 
			bkg.get_grid_y()[yidx+1], bkg.get_grid_z()[zidx+1]);
			
		// Vertex #4 @ (xidx, yidx, zidx+1)
		auto [v4x, v4y, v4z] = opts.mapc2p()(bkg.get_grid_x()[xidx], 
			bkg.get_grid_y()[yidx], bkg.get_grid_z()[zidx+1]);

		return std::array<double, 12> {v1x, v1y, v1z, v2x, v2y, v2z, 
			v3x, v3y, v3z, v4x, v4y, v4z};
	}

	// Get bounding quadrilateral in the y direction in physical coordinates
	// Must go around the edge of the grid, no cutting across the middle.
	std::array <double, 12> get_y_bound_vertices(
		const Background::Background& bkg, const Options::Options& opts,
		const int xidx, const int yidx, const int zidx)
	{
		// Vertex #1 @ (xidx, yidx, zidx)
		auto [v1x, v1y, v1z] = opts.mapc2p()(bkg.get_grid_x()[xidx], 
			bkg.get_grid_y()[yidx], bkg.get_grid_z()[zidx]);

		// Vertex #2 @ (xidx+1, yidx, zidx)
		auto [v2x, v2y, v2z] = opts.mapc2p()(bkg.get_grid_x()[xidx+1], 
			bkg.get_grid_y()[yidx], bkg.get_grid_z()[zidx]);

		// Vertex #3 @ (xidx+1, yidx, zidx+1)
		auto [v3x, v3y, v3z] = opts.mapc2p()(bkg.get_grid_x()[xidx+1], 
			bkg.get_grid_y()[yidx], bkg.get_grid_z()[zidx+1]);
			
		// Vertex #4 @ (xidx, yidx, zidx+1)
		auto [v4x, v4y, v4z] = opts.mapc2p()(bkg.get_grid_x()[xidx], 
			bkg.get_grid_y()[yidx], bkg.get_grid_z()[zidx+1]);

		return std::array<double, 12> {v1x, v1y, v1z, v2x, v2y, v2z, 
			v3x, v3y, v3z, v4x, v4y, v4z};
	}

	// Get bounding quadrilateral in the z direction in physical coordinates
	// Must go around the edge of the grid, no cutting across the middle.
	std::array <double, 12> get_z_bound_vertices(
		const Background::Background& bkg, const Options::Options& opts,
		const int xidx, const int yidx, const int zidx)
	{
		// Vertex #1 @ (xidx, yidx, zidx)
		auto [v1x, v1y, v1z] = opts.mapc2p()(bkg.get_grid_x()[xidx], 
			bkg.get_grid_y()[yidx], bkg.get_grid_z()[zidx]);

		// Vertex #2 @ (xidx+1, yidx, zidx)
		auto [v2x, v2y, v2z] = opts.mapc2p()(bkg.get_grid_x()[xidx+1], 
			bkg.get_grid_y()[yidx], bkg.get_grid_z()[zidx]);

		// Vertex #3 @ (xidx+1, yidx+1, zidx)
		auto [v3x, v3y, v3z] = opts.mapc2p()(bkg.get_grid_x()[xidx+1], 
			bkg.get_grid_y()[yidx+1], bkg.get_grid_z()[zidx]);
			
		// Vertex #4 @ (xidx, yidx+1, zidx)
		auto [v4x, v4y, v4z] = opts.mapc2p()(bkg.get_grid_x()[xidx], 
			bkg.get_grid_y()[yidx+1], bkg.get_grid_z()[zidx]);

		return std::array<double, 12> {v1x, v1y, v1z, v2x, v2y, v2z, 
			v3x, v3y, v3z, v4x, v4y, v4z};
	}

	bool check_boundary(const Background::Background& bkg, 
		Impurity& imp, const Options::Options& opts, const int tidx, 
		const int xidx, const int yidx, const int zidx, 
		const double imp_time_step, std::unique_ptr<KDTree::KDTree_t>& kdtree)
	{
		// If we've run past the time range covered by the background plasma
		// then we're done
		if (imp.get_t() > bkg.get_t_max()) return false; 

		// If we're at the edge of the computational grid in x, check for
		// crossing the boundary. Absorbing boundary condition.
		if (xidx == 0 || xidx == std::ssize(bkg.get_x()) - 1)
		{
			// Need to pass in either x+1 or x depending on condition!
			std::array <double, 12> xbv {};
			if (xidx == 0) xbv = get_x_bound_vertices(bkg, opts, xidx, yidx, 
				zidx);
			else if (xidx == std::ssize(bkg.get_x()) - 1) xbv = 
				get_x_bound_vertices(bkg, opts, xidx+1, yidx, zidx);

			// See at what fraction of the particle path we hit the boundary,
			// if we did. Breaking my formatting rules since this is a very
			// long function call (it's straightforward though).
			// Return -1 if not intersect, and the frac between 0-1 if so.
			double intersect_frac {MoellerTrumbore::get_intersect_frac(
				imp.get_prevX(), imp.get_prevY(), imp.get_prevZ(),   // p1
					imp.get_X(),     imp.get_Y(),     imp.get_Z(),   // p2 
						 xbv[0],          xbv[1],          xbv[2],   // v1
						 xbv[3],          xbv[4],          xbv[5],   // v2
						 xbv[6],          xbv[7],          xbv[8],   // v3
						 xbv[9],          xbv[10],        xbv[11])}; // v4

			// If -1, no intersect. If not, we intersected. Don't care where
			// since this is an absorbing boundary, it suffices to know we did.
			if (intersect_frac < 0)
			{
				//std::cout << "X absorbed!\n " 
				//	<< imp.get_prevX() << "    " << imp.get_X() << "\n"
				//	<< imp.get_prevY() << "    " << imp.get_Y() << "\n"
				//	<< imp.get_prevZ() << "    " << imp.get_Z() << "\n";
				return false;
			}
		}

		// If we're at the edge of the computational grid in z, check for
		// crossing the boundary. Absorbing boundary condition.
		if (zidx == 0 || zidx == std::ssize(bkg.get_z()) - 1)
		{
			std::array <double, 12> zbv {};
			if (zidx == 0) zbv = get_z_bound_vertices(bkg, opts, xidx, yidx, 
				zidx);
			else if (zidx == std::ssize(bkg.get_z()) - 1) zbv = 
				get_z_bound_vertices(bkg, opts, xidx, yidx, zidx+1);

			// See at what fraction of the particle path we hit the boundary,
			// if we did. Breaking my formatting rules since this is a very
			// long function call (it's straightforward though).
			// Return -1 if not intersect, and the frac between 0-1 if so.
			double intersect_frac {MoellerTrumbore::get_intersect_frac(
				imp.get_prevX(), imp.get_prevY(), imp.get_prevZ(),   // p1
					imp.get_X(),     imp.get_Y(),     imp.get_Z(),   // p2 
						 zbv[0],          zbv[1],          zbv[2],   // v1
						 zbv[3],          zbv[4],          zbv[5],   // v2
						 zbv[6],          zbv[7],          zbv[8],   // v3
						 zbv[9],          zbv[10],        zbv[11])}; // v4

			// If -1, no intersect. If not, we intersected. Don't care where
			// since this is an absorbing boundary, it suffices to know we did.
			if (intersect_frac < 0)
			{
				//std::cout << "Z absorbed!\n " 
				//	<< imp.get_prevX() << "    " << imp.get_X() << "\n"
				//	<< imp.get_prevY() << "    " << imp.get_Y() << "\n"
				//	<< imp.get_prevZ() << "    " << imp.get_Z() << "\n";
				return false;
			}
		}

		// If we're at the edge of the computational grid in y, check for
		// crossing the boundary. Periodic boundary condition.
		if (yidx == 0 || yidx == std::ssize(bkg.get_y()) - 1)
		{
			//std::cout << "Checking for y bound intersection..." << yidx 
			//	<< '\n';
			std::array <double, 12> ybv {};
			if (yidx == 0) ybv = get_y_bound_vertices(bkg, opts, xidx, yidx, 
				yidx);
			else if (yidx == std::ssize(bkg.get_y()) - 1) ybv = 
				get_y_bound_vertices(bkg, opts, xidx, yidx+1, zidx);

			// See at what fraction of the particle path we hit the boundary,
			// if we did. Breaking my formatting rules since this is a very
			// long function call (it's straightforward though).
			// Return -1 if not intersect, and the frac between 0-1 if so.
			double intersect_frac {MoellerTrumbore::get_intersect_frac(
				imp.get_prevX(), imp.get_prevY(), imp.get_prevZ(),   // p1
					imp.get_X(),     imp.get_Y(),     imp.get_Z(),   // p2 
						 ybv[0],          ybv[1],          ybv[2],   // v1
						 ybv[3],          ybv[4],          ybv[5],   // v2
						 ybv[6],          ybv[7],          ybv[8],   // v3
						 ybv[9],          ybv[10],        ybv[11])}; // v4

			// If -1, no intersect. If not, we intersected. Don't care where
			// since this is an absorbing boundary, it suffices to know we did.
			if (intersect_frac < 0) return true;

			// If there is an intersection, move particle back in time to put 
			// it on the boundary by subtracting the portion of the vector 
			// past the boundary. The false is set_prev. We don't want to
			// assign the previous coordinate here since we are effectively
			// taking a step backwards to do this step again in stages.
			//std::cout << " Backing up X,Y,Z:\n";
			double X_excess {(1 - intersect_frac) * (imp.get_X() 
				- imp.get_prevX())};
			double Y_excess {(1 - intersect_frac) * (imp.get_Y() 
				- imp.get_prevY())};
			double Z_excess {(1 - intersect_frac) * (imp.get_Z() 
				- imp.get_prevZ())};
			//imp.set_X(imp.get_X() - X_excess, false);
			//imp.set_Y(imp.get_Y() - Y_excess, false);
			//imp.set_Z(imp.get_Z() - Z_excess, false);

			// Additonal steps for a periodic boundary. It very difficult to
			// perfectly preserve the physical coordinates without actually
			// tracking the position in computational space (we would need
			// a mapp2c in addition to a mapc2p, but it doesn't exist). So we
			// do an approximation based on the grid nodes, since we can easily
			// say, "this node corresponds to this node on the other side of
			// "the grid". If we need to loop the particle back 
			// around the grid for periodic conditions, we first find the 
			// nearest node and move the particle there. Then we know where
			// to loop it around: if yidx = 0, then reassign yidx to 
			// grid_y.size()-1, and vice-versa. Then change the particle's
			// physical coordinates to this new node, and finish the time step
			// (potentially triggering a boundary condition again recursively).

			// Periodic boundary condition, so just wrap it around in the 
			// y index. The -1 is because of zero-indexing.
			int len_y {static_cast<int>(std::ssize(bkg.get_y()))};
			int new_yidx {};
			if (yidx == len_y - 1) new_yidx = 0;
			else if (yidx == 0) new_yidx = len_y - 1;

			// Find nearest grid node. We are checking against the four 
			// coordinates represented by the cell face we are crossing.
			double vert_dist {std::numeric_limits<double>::max()};
			int near_id {};
			for (int i {}; i < 4; i++)
			{
				// Calculate distance to each vertex, and assign its number
				// to near_id (near_id=1 is v1, etc.).
				double dist_sq {
					(imp.get_X() - ybv[i*3]) * (imp.get_X() - ybv[i*3]) +
					(imp.get_Y() - ybv[i*3+1]) * (imp.get_Y() - ybv[i*3+1]) +
					(imp.get_Z() - ybv[i*3+2]) * (imp.get_Z() - ybv[i*3+2])};
				if (dist_sq < vert_dist)
				{
					vert_dist = dist_sq;
					near_id = i + 1;
				}
			}
			
			// Reset particle position to one of these indices
			// Vertex #1 @ (xidx, yidx, zidx)
			// Vertex #2 @ (xidx+1, yidx, zidx)
			// Vertex #3 @ (xidx+1, yidx, zidx+1)
			// Vertex #4 @ (xidx, yidx, zidx+1)
			int new_xidx {};
			int new_zidx {};
			if (near_id == 1)
			{
				new_xidx = xidx;
				new_zidx = zidx;
			}
			else if (near_id == 2)
			{
				new_xidx = xidx + 1;
				new_zidx = zidx;
			}
			else if (near_id == 3)
			{
				new_xidx = xidx + 1;
				new_zidx = zidx + 1;
			}
			else if (near_id == 4)
			{
				new_xidx = xidx;
				new_zidx = zidx + 1;
			}

			// Update computational coordinates, which are actually just the
			// center coordinates of whatever cell it is in (and thus use the
			// original indices for x and z).
			imp.set_x(bkg.get_x()[xidx]);
			imp.set_y(bkg.get_y()[new_yidx]);
			imp.set_z(bkg.get_z()[zidx]);

			// Update physical coordinates at the new node on the other side 
			// of the grid.
			auto [X, Y, Z] = opts.mapc2p()(bkg.get_grid_x()[new_xidx], 
				bkg.get_grid_y()[new_yidx], bkg.get_grid_z()[new_zidx]);
			imp.set_X(X, false);
			imp.set_Y(Y, false);
			imp.set_Z(Z, false);

			// Now take the length of the step that was cutoff for going out
			// of bounds, and move the particle that distance from the node
			// it has been placed at towards the cell center. We can't just
			// continue it along the vector it was going at, because we have
			// no gaurantee it'll finish in the cell, it could end up outside
			// and then it's pretty much lost for good since we do not check if
			// the physical coordinates exist within a grid element. 

			// Unit vector pointing towards cell center, calculate components
			// then normalize to unit length
			double unit_X {bkg.get_X()(new_xidx, new_yidx, new_zidx) 
				- imp.get_X()};
			double unit_Y {bkg.get_Y()(new_xidx, new_yidx, new_zidx) 
				- imp.get_Y()};
			double unit_Z {bkg.get_Z()(new_xidx, new_yidx, new_zidx) 
				- imp.get_Z()};
			double unit_mag {std::sqrt(unit_X * unit_X + unit_Y * unit_Y 
				+ unit_Z * unit_Z)};
			unit_X = unit_X / unit_mag;
			unit_Y = unit_Y / unit_mag;
			unit_Z = unit_Z / unit_mag;

			// Magnitude of the excess distance that went past
			double excess_mag {std::sqrt(X_excess * X_excess 
				+ Y_excess * Y_excess + Z_excess * Z_excess)};

			// If excess_mag > unit_mag, we are overshooting the center
			// coordinates. This could be fine, but we don't know if we are
			// overshooting the center and blasting off the grid. Solution is
			// to decrease time step.
			if (excess_mag > unit_mag)
			{
				std::cout << "Warning! The periodic boundary algorithm has "
					<< "detected too large a time step. Consider reducing "
					<< "the timestep if possible. Debug: " << excess_mag << " " 
					<< unit_mag << '\n';
			}

			// Update particle position to move along the unit vector towards
			// the center by excess_mag distance
			imp.set_X(imp.get_X() + unit_X * excess_mag);
			imp.set_Y(imp.get_Y() + unit_Y * excess_mag);
			imp.set_Z(imp.get_Z() + unit_Z * excess_mag);
		}

		// No boundary condition met
		return true;
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
		const std::vector<int> imp_var_reduct_counts, 
		const bool imp_var_reduct_on, 
		std::unique_ptr<KDTree::KDTree_t>& kdtree,
		const Options::Options& opts)
	{
		// Timestep of impurity transport simulation. It can be a constant
		// value (set here), or set on the fly based on a reasonable criteria.
		// Initialize it with whatever the minimum time step is.
		//double imp_time_step {Input::get_opt_dbl(Input::imp_time_step_min)};
		double imp_time_step {opts.imp_time_step_min()};
		if (opts.imp_time_step_opt_int() == 0)
		{
			imp_time_step = opts.imp_time_step();
		}

		// For debugging purposes (only works with one thread)
		//static int imp_id {0};
		//imp_id++;

		bool continue_following {true};
		while (continue_following)
		{
			// Get nearest time index
			int tidx {get_nearest_index(bkg.get_times(), imp.get_t())};

			// Get x,y,z indices in computational space
			int xidx {get_nearest_cell_index(bkg.get_grid_x(), imp.get_x())};
			int yidx {get_nearest_cell_index(bkg.get_grid_y(), imp.get_y())};
			int zidx {get_nearest_cell_index(bkg.get_grid_z(), imp.get_z())};

			// Calculate Lorentz force components. First loop these are all
			// zero if particles start at rest.
			auto [fX, fY, fZ] = lorentz_forces(imp, bkg, tidx, xidx, yidx, 
				zidx);

			// Variance reduction scheme. This may change the particle's 
			// weight and create a secondary Impurity to be followed. It
			// currently is based on ionization/recombination, hence that also
			// needs to be on to work.
			if (imp_var_reduct_on && opts.imp_iz_recomb_int() > 0)
			{
				VarianceReduction::split_particle_main(imp, tidx, xidx, yidx, 
					zidx, imp_stats, imps, opts.imp_var_reduct_min_weight(), 
					imp_var_reduct_counts, bkg, oa_ioniz, oa_recomb, 
					imp_time_step);
			}

			// Calculate variable time step (if necessary)
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
			//std::cout << "Before step\n";
			//std::cout << "  id, tidx, q, t, x, y, z, dt, fX, fY, fZ, " << 
			//	"vX, vY, vZ, X, Y, Z: \n"
			//	<< " " << imp_id << ", " << tidx << ", "
			//	<< imp.get_charge() << ", "<< imp.get_t() << ", " 
			//	<< ", " << imp.get_x() << ", " << imp.get_y() 
			//	<< ", " << imp.get_z() << ", " << imp_time_step << '\n'
			//	<< "  " << fX << ", " << fY << ", " << fZ << ", " 
			//	<< imp.get_vX() << ", " << imp.get_vY() << ", " 
			//	<< imp.get_vZ() << ", "
			//	<< imp.get_X() << ", " << imp.get_Y() << ", " 
			//	<< imp.get_Z() << '\n';

			// Update statistics. Need to do this after the time step is
			// calculated (if it is), but before the particle moves into 
			// a different cell in step.
			record_stats(imp_stats, imp, bkg, tidx, xidx, yidx, zidx, 
				imp_time_step);

			// Last thing is move particle to a new location
			//step(imp, imp_time_step);
			step(imp, fX, fY, fZ, imp_time_step, bkg, kdtree);

			//std::cout << "After step\n";
			//std::cout << "  id, tidx, q, t, x, y, z, dt, fX, fY, fZ, " << 
			//	"vX, vY, vZ, X, Y, Z: \n"
			//	<< " " << imp_id << ", " << tidx << ", "
			//	<< imp.get_charge() << ", "<< imp.get_t() << ", " 
			//	<< ", " << imp.get_x() << ", " << imp.get_y() 
			//	<< ", " << imp.get_z() << ", " << imp_time_step << '\n'
			//	<< "  " << fX << ", " << fY << ", " << fZ << ", " 
			//	<< imp.get_vX() << ", " << imp.get_vY() << ", " 
			//	<< imp.get_vZ() << ", "
			//	<< imp.get_X() << ", " << imp.get_Y() << ", " 
			//	<< imp.get_Z() << '\n';

			// Check for ionization or recombination
			if (opts.imp_iz_recomb_int() > 0)
			{
				OpenADAS::ioniz_recomb(imp, bkg, oa_ioniz, oa_recomb, 
					imp_time_step, tidx, xidx, yidx, zidx, ioniz_warnings, 
					recomb_warnings);
			}

			// Check for a terminating or boundary conditions
			continue_following = check_boundary(bkg, imp, opts, tidx, xidx, 
				yidx, zidx, imp_time_step, kdtree);
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
		std::unique_ptr<KDTree::KDTree_t>& kdtree,
		const Options::Options& opts)
	{
		// Variance reduction requires ionization/recombination be on
		//if (imp_var_reduct_on && !imp_iz_recomb_on)
		if (opts.imp_var_reduct_int() > 0 && opts.imp_iz_recomb_int() == 0)
		{
			std::cout << "Warning! Variance reduction requires "
				<< "imp_ioniz_recomb be set to 'yes'. Variance reduction will"
				<< " not occur for this simulation.\n";
		}

		// Vector of counts at each frame, below which is considered a 
		// low-count region (as defined in variance_reduction.get_counts).
		std::vector<int> imp_var_reduct_counts (
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
		constexpr int prog_interval {10};
	
		// Loop through one impurity at a time, tracking it from its birth
		// time/location to the end. Dynamic scheduling likely the best here 
		// since the loop times can vary widely. Declaring bkg, oa_ioniz and
		// oa_recomb as shared is redundant, but I like being explicit. 
		int ioniz_warnings {};
		int recomb_warnings {};
		int prim_imp_count {};
		int tot_imp_count {};
		int priv_count {};
		bool imp_var_reduct_control {false};
		#pragma omp parallel for schedule(dynamic) \
			firstprivate(priv_count, imp_var_reduct_counts) \
			firstprivate(imp_var_reduct_control) \
			shared(bkg, oa_ioniz, oa_recomb) \
			reduction(+: imp_stats, ioniz_warnings, recomb_warnings) \
			reduction(+: tot_imp_count)
		for (int i = 0; i < opts.imp_num(); ++i)
		{

			// Periodically determine what number of counts qualifies for 
			// splitting a particle, if necessary. 
			if (opts.imp_var_reduct_int() > 0 && priv_count > 0)
			{
				// This variable is just an integer saying "every X particles
				// followed, update the variance reduction counts". It needs to
				// be greater than zero.
				int imp_var_reduct_mod {static_cast<int>(
					static_cast<double>(opts.imp_num()) / omp_get_num_threads() 
					* opts.imp_var_reduct_freq())};
				if (imp_var_reduct_mod == 0) imp_var_reduct_mod = 1;

				// Update variance reduction counts according to the variance
				// reduction frequency (imp_var_reduct_freq), which is a number 
				// between 0 and 1. So for example, a value of 0.10 would 
				// execute this if statement 10 times during the simulation.
				if (priv_count % imp_var_reduct_mod == 0 && priv_count > 0)
				{
					imp_var_reduct_counts = 
						VarianceReduction::get_counts(imp_stats, 1.0);

					// We use this separate boolean (initially false) to
					// control if variance reduction occurs so that the
					// simulation can run through for a bit and learn
					// what a low-count region qualifies as. 
					imp_var_reduct_control = true;

					//for (auto c : imp_var_reduct_counts)
					//{
					//	std::cout << "imp_var_reduct_counts = " << c << '\n';
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
					imp_var_reduct_counts, imp_var_reduct_control, kdtree, 
					opts);
				//std::cout << "imps.size() = " << imps.size() << '\n';
				++tot_imp_count;
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
							<< opts.imp_num() << " impurities (" 
							<< static_cast<int>(perc_complete) << "%)\n";
					}
				}
			}

			// Increment counter private to each thread
			priv_count++;
		}

		// Let user know how many, if any, warnings occured.
		print_ioniz_recomb_warn(ioniz_warnings, recomb_warnings);

		// Let user know how many secondary impurities were followed.
		std::cout << "Secondary impurities followed: " << tot_imp_count 
			- prim_imp_count << '\n';
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

		// Build KD-tree, which is used for nearest-neighbor searches between
		// particle X,Y,Z position and a known X,Y,Z position of the background
		// plasma.
		std::cout << "Building KDTree...\n";
		std::unique_ptr<KDTree::KDTree_t> kdtree {
			KDTree::build_tree(bkg.get_X(), bkg.get_Y(), bkg.get_Z())};

		// Execute main particle following loop.
		main_loop(bkg, imp_stats, oa_ioniz, oa_recomb, kdtree, opts);
	
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
