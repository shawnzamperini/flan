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

		// We need to start the impurity at a grid node, otherwise the 
		// code may autolocate it somewhere further away. This is actually
		// super subtle yet extremely important. I lost years off my life
		// figuring this out.
		int xidx {get_nearest_cell_index(bkg.get_grid_x(), start_x)};
		//std::cout << "start xidx: " << xidx << '\n';
		return bkg.get_grid_x()[xidx];

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
		//std::cout << "start yidx: " << yidx << '\n';
		return bkg.get_grid_y()[yidx];
	}

	double get_birth_z(const Background::Background& bkg,
		const Options::Options& opts)
	{
		double start_z {opts.imp_zstart_val()};
		int zidx {get_nearest_cell_index(bkg.get_grid_z(), start_z)};
		//std::cout << "start zidx: " << zidx << '\n';
		return bkg.get_grid_z()[zidx];
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
	
		//double R {std::sqrt(X_imp*X_imp + Y_imp*Y_imp)};
		//std::cout << "Primary created at: " << X_imp << ", " << Y_imp << 
		//	", " << Z_imp << "  (" << x_imp << ", " << y_imp << ", "
		//	<< z_imp << ")  R = " << R << "\n";

		// It is prudent to do a sanity check here that the cell checking
		// algorithm registers that the particle is in the cell it starts in.
		// If this error trips, something is wrong with the algorithm.
		int xidx {get_nearest_cell_index(bkg.get_grid_x(), x_imp)};
		int yidx {get_nearest_cell_index(bkg.get_grid_y(), y_imp)};
		int zidx {get_nearest_cell_index(bkg.get_grid_z(), z_imp)};
		bool in_cell {check_in_cell(bkg, X_imp, Y_imp, Z_imp, xidx, yidx, 
			zidx, false)};
		if (!in_cell)
		{
			std::cerr << "Error! Impurity failed sanity check on if the "
				<< "cell checking algorithm registers the intial starting cell"
				<< " correctly. This is a programming error and must be fixed."
				<< '\n' << "  (xidx,yidx,zidx) = (" << xidx << ", " << yidx 
				<< ", " << zidx << ")  (X,Y,Z) = (" << X_imp << ", " << Y_imp
				<< ", " << Z_imp << ")\n";
		}

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

	// Returns true if particle intersects x boundary
	bool check_boundary_x(const Impurity& imp, 
		const Background::Background& bkg)
	{

		// Need to check for crossing at any of the bounding elements at the
		// x boundary, not just the cell the particle is in. This is because
		// the particle could feasibly hit a neighboring x-boundary. To be
		// robust, check all possibilities.

		// Array to hold vertices of bounding quadrilateral and boolean for
		// determining if particle crossed boundary
		std::array <double, 12> bv {};
		bool intersect {false};

		for (int j {}; j < bkg.get_X().get_dim2(); j++)
		{
			for (int k {}; k < bkg.get_X().get_dim3(); k++)
			{

				// Check lower bound at xidx=0
				bv = get_x_bound_vertices(bkg, 0, j, k);
				intersect = MoellerTrumbore::check_intersect(
					imp.get_prevX(),  imp.get_prevY(),  imp.get_prevZ(), // p1
					imp.get_X(),      imp.get_Y(),      imp.get_Z(),     // p2 
					bv[0],            bv[1],            bv[2],           // v1
					bv[3],            bv[4],            bv[5],           // v2
					bv[6],            bv[7],            bv[8],           // v3
					bv[9],            bv[10],           bv[11]);         // v4

				// If intersect, then we crossed an absorbing boundary so 
				// return false
				if (intersect)
				{
					std::cout << "Absorbed x (lower)\n";
					return true;
				}

				// Check upper bound.
				bv = get_x_bound_vertices(bkg, bkg.get_X().get_dim1(), j, k);
				intersect = MoellerTrumbore::check_intersect(
					imp.get_prevX(),  imp.get_prevY(),  imp.get_prevZ(), // p1
					imp.get_X(),      imp.get_Y(),      imp.get_Z(),     // p2 
					bv[0],            bv[1],            bv[2],           // v1
					bv[3],            bv[4],            bv[5],           // v2
					bv[6],            bv[7],            bv[8],           // v3
					bv[9],            bv[10],           bv[11]);         // v4

				// If intersect, then we crossed an absorbing boundary so 
				// return false
				if (intersect)
				{
					std::cout << "Absorbed x (upper)\n";
					return true;
				}
			}
		}

		// No intersection found
		return false;
	}


	void check_boundary_y(Impurity& imp, const Background::Background& bkg)
	{
		// Temporary y boundary check that works for simple helical
		if (imp.get_Z() < bkg.get_grid_y()[0])
		{
			imp.set_Z(bkg.get_grid_y().back() + (imp.get_Z() 
				- bkg.get_grid_y()[0]));
			//std::cout << "Relected (lower)\n";
		}
		else if (imp.get_Z() > bkg.get_grid_y().back())
		{
			imp.set_Z(bkg.get_grid_y()[0] + (imp.get_Z() 
				- bkg.get_grid_y().back()));
			//std::cout << "Relected (upper)\n";
		}

	}

	// Returns true if particle intersects z boundary
	bool check_boundary_z(const Impurity& imp, 
		const Background::Background& bkg)
	{

		// Need to check for crossing at any of the bounding elements at the
		// z boundary, not just the cell the particle is in. This is because
		// the particle could feasibly hit a neighboring z-boundary. To be
		// robust, check all possibilities.

		// Array to hold vertices of bounding quadrilateral and boolean for
		// determining if particle crossed boundary
		std::array <double, 12> bv {};
		bool intersect {false};

		for (int i {}; i < bkg.get_X().get_dim1(); i++)
		{
			for (int j {}; j < bkg.get_X().get_dim2(); j++)
			{

				// Check lower bound at xidx=0
				bv = get_z_bound_vertices(bkg, i, j, 0);
				intersect = MoellerTrumbore::check_intersect(
					imp.get_prevX(),  imp.get_prevY(),  imp.get_prevZ(), // p1
					imp.get_X(),      imp.get_Y(),      imp.get_Z(),     // p2 
					bv[0],            bv[1],            bv[2],           // v1
					bv[3],            bv[4],            bv[5],           // v2
					bv[6],            bv[7],            bv[8],           // v3
					bv[9],            bv[10],           bv[11]);         // v4

				// If intersect, then we crossed an absorbing boundary so 
				// return false
				if (intersect) return true;

				// Check upper bound.
				bv = get_z_bound_vertices(bkg, i, j, bkg.get_X().get_dim3());
				intersect = MoellerTrumbore::check_intersect(
					imp.get_prevX(),  imp.get_prevY(),  imp.get_prevZ(), // p1
					imp.get_X(),      imp.get_Y(),      imp.get_Z(),     // p2 
					bv[0],            bv[1],            bv[2],           // v1
					bv[3],            bv[4],            bv[5],           // v2
					bv[6],            bv[7],            bv[8],           // v3
					bv[9],            bv[10],           bv[11]);         // v4

				// If intersect, then we crossed an absorbing boundary so 
				// return false
				if (intersect) return true;
			}
		}

		// No intersection found
		return false;
	}

	bool step(Impurity& imp, const double fX, const double fY, const double fZ, 
		const double imp_time_step, const Background::Background& bkg,
		std::unique_ptr<KDTree::KDTree_t>& kdtree, 
		const Options::Options& opts, const int tidx, int& xidx, int& yidx, 
		int& zidx, bool& continue_following)
	{
		// Change in velocity over time step (this is just F = m * dv/dt)
		double dvX {fX * imp_time_step / imp.get_mass()};
		double dvY {fY * imp_time_step / imp.get_mass()};
		double dvZ {fZ * imp_time_step / imp.get_mass()};

		// Update particle velocity
		imp.set_vX(imp.get_vX() + dvX);
		imp.set_vY(imp.get_vY() + dvY);
		imp.set_vZ(imp.get_vZ() + dvZ);

		// Update particle time and position in physical space, one
		// coordinate a a time. 

		// If we've run past the time range covered by the background plasma
		// then we're done
		imp.set_t(imp.get_t() + imp_time_step);
		if (imp.get_t() > bkg.get_t_max())
		{
			return false;
		}

		imp.set_X(imp.get_X() + imp.get_vX() * imp_time_step);
		imp.set_Y(imp.get_Y() + imp.get_vY() * imp_time_step);
		imp.set_Z(imp.get_Z() + imp.get_vZ() * imp_time_step);

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
	//std::array <double, 12> get_x_bound_vertices(
	//	const Background::Background& bkg, const Options::Options& opts,
	//	const int xidx, const int yidx, const int zidx)
	std::array <double, 12> get_x_bound_vertices(
		const Background::Background& bkg,
		const int xidx, const int yidx, const int zidx)
	{
		//std::cout << "replace with grid_X!\n";
		// Vertex #1 @ (xidx, yidx, zidx)
		//auto [v1x, v1y, v1z] = opts.mapc2p()(bkg.get_grid_x()[xidx], 
		//	bkg.get_grid_y()[yidx], bkg.get_grid_z()[zidx]);
		double v1x {bkg.get_grid_X()(xidx, yidx, zidx)};
		double v1y {bkg.get_grid_Y()(xidx, yidx, zidx)};
		double v1z {bkg.get_grid_Z()(xidx, yidx, zidx)};

		// Vertex #2 @ (xidx, yidx+1, zidx)
		//auto [v2x, v2y, v2z] = opts.mapc2p()(bkg.get_grid_x()[xidx], 
		//	bkg.get_grid_y()[yidx+1], bkg.get_grid_z()[zidx]);
		double v2x {bkg.get_grid_X()(xidx, yidx+1, zidx)};
		double v2y {bkg.get_grid_Y()(xidx, yidx+1, zidx)};
		double v2z {bkg.get_grid_Z()(xidx, yidx+1, zidx)};

		// Vertex #3 @ (xidx, yidx+1, zidx+1)
		//auto [v3x, v3y, v3z] = opts.mapc2p()(bkg.get_grid_x()[xidx], 
		//	bkg.get_grid_y()[yidx+1], bkg.get_grid_z()[zidx+1]);
		double v3x {bkg.get_grid_X()(xidx, yidx+1, zidx+1)};
		double v3y {bkg.get_grid_Y()(xidx, yidx+1, zidx+1)};
		double v3z {bkg.get_grid_Z()(xidx, yidx+1, zidx+1)};
			
		// Vertex #4 @ (xidx, yidx, zidx+1)
		//auto [v4x, v4y, v4z] = opts.mapc2p()(bkg.get_grid_x()[xidx], 
		//	bkg.get_grid_y()[yidx], bkg.get_grid_z()[zidx+1]);
		double v4x {bkg.get_grid_X()(xidx, yidx, zidx+1)};
		double v4y {bkg.get_grid_Y()(xidx, yidx, zidx+1)};
		double v4z {bkg.get_grid_Z()(xidx, yidx, zidx+1)};

		return std::array<double, 12> {v1x, v1y, v1z, v2x, v2y, v2z, 
			v3x, v3y, v3z, v4x, v4y, v4z};
	}

	// Get bounding quadrilateral in the y direction in physical coordinates
	// Must go around the edge of the grid, no cutting across the middle.
	std::array <double, 12> get_y_bound_vertices(
		const Background::Background& bkg,
		const int xidx, const int yidx, const int zidx)
	{
		// Vertex #1 @ (xidx, yidx, zidx)
		//auto [v1x, v1y, v1z] = opts.mapc2p()(bkg.get_grid_x()[xidx], 
		//	bkg.get_grid_y()[yidx], bkg.get_grid_z()[zidx]);
		double v1x {bkg.get_grid_X()(xidx, yidx, zidx)};
		double v1y {bkg.get_grid_Y()(xidx, yidx, zidx)};
		double v1z {bkg.get_grid_Z()(xidx, yidx, zidx)};

		// Vertex #2 @ (xidx+1, yidx, zidx)
		//auto [v2x, v2y, v2z] = opts.mapc2p()(bkg.get_grid_x()[xidx+1], 
		//	bkg.get_grid_y()[yidx], bkg.get_grid_z()[zidx]);
		double v2x {bkg.get_grid_X()(xidx+1, yidx, zidx)};
		double v2y {bkg.get_grid_Y()(xidx+1, yidx, zidx)};
		double v2z {bkg.get_grid_Z()(xidx+1, yidx, zidx)};

		// Vertex #3 @ (xidx+1, yidx, zidx+1)
		//auto [v3x, v3y, v3z] = opts.mapc2p()(bkg.get_grid_x()[xidx+1], 
		//	bkg.get_grid_y()[yidx], bkg.get_grid_z()[zidx+1]);
		double v3x {bkg.get_grid_X()(xidx+1, yidx, zidx+1)};
		double v3y {bkg.get_grid_Y()(xidx+1, yidx, zidx+1)};
		double v3z {bkg.get_grid_Z()(xidx+1, yidx, zidx+1)};
			
		// Vertex #4 @ (xidx, yidx, zidx+1)
		//auto [v4x, v4y, v4z] = opts.mapc2p()(bkg.get_grid_x()[xidx], 
		//	bkg.get_grid_y()[yidx], bkg.get_grid_z()[zidx+1]);
		double v4x {bkg.get_grid_X()(xidx, yidx, zidx+1)};
		double v4y {bkg.get_grid_Y()(xidx, yidx, zidx+1)};
		double v4z {bkg.get_grid_Z()(xidx, yidx, zidx+1)};

		return std::array<double, 12> {v1x, v1y, v1z, v2x, v2y, v2z, 
			v3x, v3y, v3z, v4x, v4y, v4z};
	}

	// Get bounding quadrilateral in the z direction in physical coordinates
	// Must go around the edge of the grid, no cutting across the middle.
	std::array <double, 12> get_z_bound_vertices(
		const Background::Background& bkg,
		const int xidx, const int yidx, const int zidx)
	{
		// Vertex #1 @ (xidx, yidx, zidx)
		//auto [v1x, v1y, v1z] = opts.mapc2p()(bkg.get_grid_x()[xidx], 
		//	bkg.get_grid_y()[yidx], bkg.get_grid_z()[zidx]);
		double v1x {bkg.get_grid_X()(xidx, yidx, zidx)};
		double v1y {bkg.get_grid_Y()(xidx, yidx, zidx)};
		double v1z {bkg.get_grid_Z()(xidx, yidx, zidx)};

		// Vertex #2 @ (xidx+1, yidx, zidx)
		//auto [v2x, v2y, v2z] = opts.mapc2p()(bkg.get_grid_x()[xidx+1], 
		//	bkg.get_grid_y()[yidx], bkg.get_grid_z()[zidx]);
		double v2x {bkg.get_grid_X()(xidx+1, yidx, zidx)};
		double v2y {bkg.get_grid_Y()(xidx+1, yidx, zidx)};
		double v2z {bkg.get_grid_Z()(xidx+1, yidx, zidx)};

		// Vertex #3 @ (xidx+1, yidx+1, zidx)
		//auto [v3x, v3y, v3z] = opts.mapc2p()(bkg.get_grid_x()[xidx+1], 
		//	bkg.get_grid_y()[yidx+1], bkg.get_grid_z()[zidx]);
		double v3x {bkg.get_grid_X()(xidx+1, yidx+1, zidx)};
		double v3y {bkg.get_grid_Y()(xidx+1, yidx+1, zidx)};
		double v3z {bkg.get_grid_Z()(xidx+1, yidx+1, zidx)};
			
		// Vertex #4 @ (xidx, yidx+1, zidx)
		//auto [v4x, v4y, v4z] = opts.mapc2p()(bkg.get_grid_x()[xidx], 
		//	bkg.get_grid_y()[yidx+1], bkg.get_grid_z()[zidx]);
		double v4x {bkg.get_grid_X()(xidx, yidx+1, zidx)};
		double v4y {bkg.get_grid_Y()(xidx, yidx+1, zidx)};
		double v4z {bkg.get_grid_Z()(xidx, yidx+1, zidx)};

		return std::array<double, 12> {v1x, v1y, v1z, v2x, v2y, v2z, 
			v3x, v3y, v3z, v4x, v4y, v4z};
	}

	bool check_in_cell(const Background::Background& bkg, 
		const double X, const double Y, const double Z, const int xidx, 
		const int yidx, const int zidx, const bool debug)
	{
		// Some variables used below
		int num_intersects {0};
		bool intersect {false};

		// The bounds vertices (bv)
		std::array <double, 12> bv {};

		// The way this algorithm works is by drawing a ray from X,Y,Z
		// to the origin. Going to the origin is arbitrary,
		// it just needs to sart at X, Y, Z. Then we count how many times
		// the ray intersects the surfaces of the cell at xidx, yidx, zidx.
		// If there is an odd number of intersects, then X, Y, Z
		// is inside the cell. If there is an even amount (including 0), then
		// X, Y, Z is outside the cell. So to test this, we check for 
		// intersections on each of the 6 faces of the cell.
		// I had a tough time with this - p1 needs to be X,Y,Z and p2 the
		// the origin. It will not work if you have these swapped!

		// Check lower x bounds
		bv = get_x_bound_vertices(bkg, xidx, yidx, zidx);
		intersect = MoellerTrumbore::check_intersect(
	            X,              Y,              Z,   // p1 
              0.0,            0.0,              Z,   // p2
			bv[0],          bv[1],          bv[2],   // v1
			bv[3],          bv[4],          bv[5],   // v2
			bv[6],          bv[7],          bv[8],   // v3
			bv[9],         bv[10],         bv[11], debug);  // v4
		if (intersect)
		{
			num_intersects++;
		}

		// Check upper x bounds
		bv = get_x_bound_vertices(bkg, xidx+1, yidx, zidx);
		intersect = MoellerTrumbore::check_intersect(
	            X,              Y,              Z,   // p1 
              0.0,            0.0,              Z,   // p2
			bv[0],          bv[1],          bv[2],   // v1
			bv[3],          bv[4],          bv[5],   // v2
			bv[6],          bv[7],          bv[8],   // v3
			bv[9],         bv[10],         bv[11], debug);  // v4
		if (intersect)
		{
			num_intersects++;
		}

		// Check lower y bounds
		bv = get_y_bound_vertices(bkg, xidx, yidx, zidx);
		intersect = MoellerTrumbore::check_intersect(
	            X,              Y,              Z,   // p1 
              0.0,            0.0,              Z,   // p2
			bv[0],          bv[1],          bv[2],   // v1
			bv[3],          bv[4],          bv[5],   // v2
			bv[6],          bv[7],          bv[8],   // v3
			bv[9],         bv[10],         bv[11], debug);  // v4
		if (intersect)
		{
			num_intersects++;
		}

		// Check upper y bounds
		bv = get_y_bound_vertices(bkg, xidx, yidx+1, zidx);
		intersect = MoellerTrumbore::check_intersect(
	            X,              Y,              Z,   // p1 
              0.0,            0.0,              Z,   // p2
			bv[0],          bv[1],          bv[2],   // v1
			bv[3],          bv[4],          bv[5],   // v2
			bv[6],          bv[7],          bv[8],   // v3
			bv[9],         bv[10],         bv[11], debug);  // v4
		if (intersect)
		{
			num_intersects++;
		}

		// Check lower z bounds
		bv = get_z_bound_vertices(bkg, xidx, yidx, zidx);
		intersect = MoellerTrumbore::check_intersect(
	            X,              Y,              Z,   // p1 
              0.0,            0.0,              Z,   // p2
			bv[0],          bv[1],          bv[2],   // v1
			bv[3],          bv[4],          bv[5],   // v2
			bv[6],          bv[7],          bv[8],   // v3
			bv[9],         bv[10],         bv[11], debug);  // v4
		if (intersect)
		{
			num_intersects++;
		}

		// Check upper z bounds
		bv = get_z_bound_vertices(bkg, xidx, yidx, zidx+1);
		intersect = MoellerTrumbore::check_intersect(
	            X,              Y,              Z,   // p1 
              0.0,            0.0,              Z,   // p2
			bv[0],          bv[1],          bv[2],   // v1
			bv[3],          bv[4],          bv[5],   // v2
			bv[6],          bv[7],          bv[8],   // v3
			bv[9],         bv[10],         bv[11], debug);  // v4
		if (intersect)
		{
			num_intersects++;
		}

		// Even intersections = not in cell
		//std::cout << "num_intersects = " << num_intersects << '\n';
		if (num_intersects % 2 == 0) return false;
		else
		{
			/*
			std::cout << "In cell! num_intersects = " << num_intersects << '\n';
			std::cout << X << ", " << Y << ", " << Z << '\n';	
			std::cout << "---" << xidx << " " << yidx << " " << zidx << "---\n";	
			bv = get_x_bound_vertices(bkg, xidx, yidx, zidx);
			for (auto b : bv) std::cout << b << " ";
			std::cout << '\n';
			std::cout << "  R = " 
				<< std::sqrt(bv[0]*bv[0] + bv[1]*bv[1]) << " " 
				<< std::sqrt(bv[3]*bv[3] + bv[4]*bv[4]) << " " 
				<< std::sqrt(bv[6]*bv[6] + bv[7]*bv[7]) << " " 
				<< std::sqrt(bv[9]*bv[9] + bv[10]*bv[10]) << '\n'; 
			intersect = MoellerTrumbore::check_intersect(
					X,              Y,              Z,   // p1 
				  0.0,            0.0,              Z,   // p2
				bv[0],          bv[1],          bv[2],   // v1
				bv[3],          bv[4],          bv[5],   // v2
				bv[6],          bv[7],          bv[8],   // v3
				bv[9],         bv[10],         bv[11], true);  // v4

			bv = get_x_bound_vertices(bkg, xidx+1, yidx, zidx);
			for (auto b : bv) std::cout << b << " ";
			std::cout << '\n';
			std::cout << "  R = " 
				<< std::sqrt(bv[0]*bv[0] + bv[1]*bv[1]) << " " 
				<< std::sqrt(bv[3]*bv[3] + bv[4]*bv[4]) << " " 
				<< std::sqrt(bv[6]*bv[6] + bv[7]*bv[7]) << " " 
				<< std::sqrt(bv[9]*bv[9] + bv[10]*bv[10]) << '\n'; 
			intersect = MoellerTrumbore::check_intersect(
					X,              Y,              Z,   // p1 
				  0.0,            0.0,              Z,   // p2
				bv[0],          bv[1],          bv[2],   // v1
				bv[3],          bv[4],          bv[5],   // v2
				bv[6],          bv[7],          bv[8],   // v3
				bv[9],         bv[10],         bv[11], true);  // v4

			bv = get_y_bound_vertices(bkg, xidx, yidx, zidx);
			for (auto b : bv) std::cout << b << " ";
			std::cout << '\n';
			std::cout << "  R = " 
				<< std::sqrt(bv[0]*bv[0] + bv[1]*bv[1]) << " " 
				<< std::sqrt(bv[3]*bv[3] + bv[4]*bv[4]) << " " 
				<< std::sqrt(bv[6]*bv[6] + bv[7]*bv[7]) << " " 
				<< std::sqrt(bv[9]*bv[9] + bv[10]*bv[10]) << '\n'; 
			intersect = MoellerTrumbore::check_intersect(
					X,              Y,              Z,   // p1 
				  0.0,            0.0,              Z,   // p2
				bv[0],          bv[1],          bv[2],   // v1
				bv[3],          bv[4],          bv[5],   // v2
				bv[6],          bv[7],          bv[8],   // v3
				bv[9],         bv[10],         bv[11], true);  // v4

			bv = get_y_bound_vertices(bkg, xidx, yidx+1, zidx);
			for (auto b : bv) std::cout << b << " ";
			std::cout << '\n';
			std::cout << "  R = " 
				<< std::sqrt(bv[0]*bv[0] + bv[1]*bv[1]) << " " 
				<< std::sqrt(bv[3]*bv[3] + bv[4]*bv[4]) << " " 
				<< std::sqrt(bv[6]*bv[6] + bv[7]*bv[7]) << " " 
				<< std::sqrt(bv[9]*bv[9] + bv[10]*bv[10]) << '\n'; 
			intersect = MoellerTrumbore::check_intersect(
					X,              Y,              Z,   // p1 
				  0.0,            0.0,              Z,   // p2
				bv[0],          bv[1],          bv[2],   // v1
				bv[3],          bv[4],          bv[5],   // v2
				bv[6],          bv[7],          bv[8],   // v3
				bv[9],         bv[10],         bv[11], true);  // v4

			bv = get_z_bound_vertices(bkg, xidx, yidx, zidx);
			for (auto b : bv) std::cout << b << " ";
			std::cout << '\n';
			std::cout << "  R = " 
				<< std::sqrt(bv[0]*bv[0] + bv[1]*bv[1]) << " " 
				<< std::sqrt(bv[3]*bv[3] + bv[4]*bv[4]) << " " 
				<< std::sqrt(bv[6]*bv[6] + bv[7]*bv[7]) << " " 
				<< std::sqrt(bv[9]*bv[9] + bv[10]*bv[10]) << '\n'; 
			intersect = MoellerTrumbore::check_intersect(
					X,              Y,              Z,   // p1 
				  0.0,            0.0,              Z,   // p2
				bv[0],          bv[1],          bv[2],   // v1
				bv[3],          bv[4],          bv[5],   // v2
				bv[6],          bv[7],          bv[8],   // v3
				bv[9],         bv[10],         bv[11], true);  // v4

			bv = get_z_bound_vertices(bkg, xidx, yidx, zidx+1);
			for (auto b : bv) std::cout << b << " ";
			std::cout << '\n';
			std::cout << "  R = " 
				<< std::sqrt(bv[0]*bv[0] + bv[1]*bv[1]) << " " 
				<< std::sqrt(bv[3]*bv[3] + bv[4]*bv[4]) << " " 
				<< std::sqrt(bv[6]*bv[6] + bv[7]*bv[7]) << " " 
				<< std::sqrt(bv[9]*bv[9] + bv[10]*bv[10]) << '\n'; 
			intersect = MoellerTrumbore::check_intersect(
					X,              Y,              Z,   // p1 
				  0.0,            0.0,              Z,   // p2
				bv[0],          bv[1],          bv[2],   // v1
				bv[3],          bv[4],          bv[5],   // v2
				bv[6],          bv[7],          bv[8],   // v3
				bv[9],         bv[10],         bv[11], true);  // v4
			std::cout << "-----\n";	
			*/
			return true;
		}
		return false;
	}
	bool find_containing_cell(Impurity& imp, 
		const Background::Background& bkg, 
		int& xidx, int& yidx, int& zidx, 
		std::unique_ptr<KDTree::KDTree_t>& kdtree)
	{

		// Instead of trying to limit the search by finding what the nearest 
		// node is, instead brute force it the first go through by checking
		// every cell. Then we can pass in the previous cell idx and search
		// outwards from that one each time to decrease time until finding.
		bool in_cell {false};

		// See if the particle is in the last known cell, or any of the ones
		// nearby since that's most likely. 
		for (int i {xidx-1}; i < xidx+2; i++)
		{
		for (int j {yidx-1}; j < yidx+2; j++)
		{
		for (int k {zidx-1}; k < zidx+2; k++)
		{
			// Additional checks to make sure we don't accidentally index 
			// outside the allowed range.
			if (i > 0 && i < bkg.get_X().get_dim1() &&
				j > 0 && j < bkg.get_X().get_dim2() &&
				k > 0 && k < bkg.get_X().get_dim3())
			{
				//std::cout << imp.get_X() << ", " << imp.get_Y() << ", " 
				//	<< imp.get_Z() << "  checking  " << i << ", " << j 
				//	<< ", " << k << '\n';
				in_cell = check_in_cell(bkg, imp.get_X(), imp.get_Y(), 
					imp.get_Z(), i, j, k);

				// If in cell, we can assign the candidate idxs to the in/out
				// parameters and return true.
				if (in_cell)
				{
					//std::cout << "In cell: " << i << ", " 
					//	<< j << ", " << k << "   " << bkg.get_x()[i] 
					//	<< ", " << bkg.get_y()[j] << ", " << bkg.get_z()[k]
					//	<< "\n";
					xidx = i;
					yidx = j;
					zidx = k;
					imp.set_x(bkg.get_x()[xidx]);
					imp.set_y(bkg.get_y()[yidx]);
					imp.set_z(bkg.get_z()[zidx]);
					return true;
				}
			}
		}
		}
		}
	
		// If you get here, the particle is more than one cell away. Not ideal,
		// and may indicate poor simulation setup but it could still happen in
		// good simulations. For instance, if the particle crossed a periodic
		// boundary.
		for (int i {}; i < bkg.get_X().get_dim1(); i++)
		{
		for (int j {}; j < bkg.get_X().get_dim2(); j++)
		{
		for (int k {}; k < bkg.get_X().get_dim3(); k++)
		{
			in_cell = check_in_cell(bkg, imp.get_X(), imp.get_Y(), imp.get_Z(),
				i, j, k);

			// If in cell, we can assign the candidate idxs to the in/out
			// parameters and return true.
			if (in_cell)
			{
				//std::cout << "In cell: " << i << ", " 
				//	<< j << ", " << k << '\n';
				xidx = i;
				yidx = j;
				zidx = k;
				imp.set_x(bkg.get_x()[xidx]);
				imp.set_y(bkg.get_y()[yidx]);
				imp.set_z(bkg.get_z()[zidx]);
				return true;
			}
		}
		}
		}
	
		return false;
	}

	bool check_boundary(const Background::Background& bkg, 
		Impurity& imp, const Options::Options& opts, const int tidx, 
		const int xidx, const int yidx, const int zidx, 
		const double imp_time_step, std::unique_ptr<KDTree::KDTree_t>& kdtree)
	{

		// Check for x boundary cross at edge of grid (absorbing)
		bool intersect_x {false};
		intersect_x = check_boundary_x(imp, bkg);

		// Check for z boundary cross at edge of grid (absorbing)
		bool intersect_z {false};
		intersect_z = check_boundary_z(imp, bkg);

		// Check for y boundary cross at edge of grid (periodic)
		check_boundary_y(imp, bkg);

		// Check for unrealistic condition indicating the algorithm isn't
		// working correctly.
		if (intersect_x && intersect_z)
		{
			std::cerr << "Error! A particle has somehow crossed multiple "
				<< "computational boundaries at the same time. This is a "
				<< "probably a programming error.\n";
		}

		// One of the absorbing boundaries has been crossed, so return false.
		if (intersect_x || intersect_z) return false;
			
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

		bool continue_following {true};

		// Get nearest time index
		int tidx {get_nearest_index(bkg.get_times(), imp.get_t())};

		// Get x,y,z indices in computational space
		int xidx {get_nearest_cell_index(bkg.get_grid_x(), imp.get_x())};
		int yidx {get_nearest_cell_index(bkg.get_grid_y(), imp.get_y())};
		int zidx {get_nearest_cell_index(bkg.get_grid_z(), imp.get_z())};
		//int xidx {};
		//int yidx {};
		//int zidx {};
		//std::cout << "Before: " << xidx << ", " << yidx << ", " << zidx << '\n';
		bool in_cell {find_containing_cell(imp, bkg, xidx, yidx, zidx, kdtree)};
		//std::cout << "After: " << xidx << ", " << yidx << ", " << zidx << '\n';

		// If not in a cell there's an issue in where the particle starts. If
		// it is in a cell, record starting position.
		if (!in_cell)
		{
			std::cerr << "Error! Particle starting outside of grid. Not "
				<< "following\n";
			continue_following = false;
		}
		else
		{
			record_stats(imp_stats, imp, bkg, tidx, xidx, yidx, zidx, 
				imp_time_step);
		}

		// For debugging purposes (only works with one thread)
		static int imp_id {0};
		imp_id++;
		bool debug {false};

		while (continue_following)
		{
			// Get nearest time index
			int tidx {get_nearest_index(bkg.get_times(), imp.get_t())};

			// Get x,y,z indices in computational space
			//int xidx {get_nearest_cell_index(bkg.get_grid_x(), imp.get_x())};
			//int yidx {get_nearest_cell_index(bkg.get_grid_y(), imp.get_y())};
			//int zidx {get_nearest_cell_index(bkg.get_grid_z(), imp.get_z())};

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
			if (imp_var_reduct_on && opts.imp_iz_recomb_int() > 0)
			{
				VarianceReduction::split_particle_main(imp, tidx, xidx, yidx, 
					zidx, imp_stats, imps, opts.imp_var_reduct_min_weight(), 
					imp_var_reduct_counts, bkg, oa_ioniz, oa_recomb, 
					imp_time_step);
			}

			// Calculate variable time step (if necessary). NOT TRUSTWORTHY
			// YET AFTER UPGRADING PARTICLE FOLLOWING.
			//if (opts.imp_time_step_opt_int() == 1) imp_time_step = 
			//	get_var_time_step(imp, bkg, tidx, xidx, yidx, zidx, fX, 
			//	fY, fZ, opts);

			// Check for a collision
			if (opts.imp_collisions_int() > 0)
			{
				collision(imp, bkg, imp_time_step, tidx, xidx, yidx, zidx, 
					opts);
			}

			// Debugging
			if (debug)
			{
				std::cout << "Before step\n";
				std::cout << "  id, tidx, q, t, x, y, z, dt, fX, fY, fZ, " << 
					"vX, vY, vZ, X, Y, Z: \n"
					<< " " << imp_id << ", " << tidx << ", "
					<< imp.get_charge() << ", "<< imp.get_t() << ", " 
					<< ", " << imp.get_x() << ", " << imp.get_y() 
					<< ", " << imp.get_z() << ", " << imp_time_step << '\n'
					<< "  " << fX << ", " << fY << ", " << fZ << ", " 
					<< imp.get_vX() << ", " << imp.get_vY() << ", " 
					<< imp.get_vZ() << ", "
					<< imp.get_X() << ", " << imp.get_Y() << ", " 
					<< imp.get_Z() << '\n';
				std::cout << tidx << " " << xidx << " " << yidx << " " 
					<< zidx << '\n';
			}

			// Perform particle step, update time and tidx. The final location 
			// of the particle could end up out of the grid. That's what the 
			// following code does before recording stats.
			continue_following = step(imp, fX, fY, fZ, imp_time_step, bkg, 
				kdtree, opts, tidx, xidx, yidx, zidx, continue_following);

			if (debug)
			{
				std::cout << "After step\n";
				std::cout << "  id, tidx, q, t, x, y, z, dt, fX, fY, fZ, " << 
					"vX, vY, vZ, X, Y, Z: \n"
					<< " " << imp_id << ", " << tidx << ", "
					<< imp.get_charge() << ", "<< imp.get_t() << ", " 
					<< ", " << imp.get_x() << ", " << imp.get_y() 
					<< ", " << imp.get_z() << ", " << imp_time_step << '\n'
					<< "  " << fX << ", " << fY << ", " << fZ << ", " 
					<< imp.get_vX() << ", " << imp.get_vY() << ", " 
					<< imp.get_vZ() << ", "
					<< imp.get_X() << ", " << imp.get_Y() << ", " 
					<< imp.get_Z() << '\n';
				std::cout << tidx << " " << xidx << " " << yidx << " " 
					<< zidx << '\n';
			}

			// Returns true if a cell was found, updating xidx, yidx, zidx 
			// internally only if a containing cell was found (i.e., the 
			// particle is in a grid cell).	
			if (debug)
			{
				//std::cout << "Finding containing cell...\n";
			}
			in_cell = find_containing_cell(imp, bkg, xidx, yidx, zidx, kdtree);
			
			// Call boundary checking algorithms if cell not found, means the 
			// particle must be outside the grid and thus crossed a boundary. 
			// The indices will be the particle's previous indices before the 
			// step since step will not have updated them.
			if (debug)
			{
				//std::cout << "Checking boundary...\n";
			}
			if (!in_cell) continue_following = check_boundary(bkg, imp, opts,
				tidx, xidx, yidx, zidx, imp_time_step, kdtree);

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
			shared(bkg, oa_ioniz, oa_recomb, kdtree) \
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
		// particle X,Y,Z position and a known X,Y,Z position of a grid node.
		std::cout << "Building KDTree...\n";
		//std::unique_ptr<KDTree::KDTree_t> kdtree {
		//	KDTree::build_tree(bkg.get_grid_X(), bkg.get_grid_Y(), 
		//	bkg.get_grid_Z())};
		std::unique_ptr<KDTree::KDTree_t> kdtree {
			KDTree::build_tree(bkg.get_X(), bkg.get_Y(), 
			bkg.get_Z())};

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
