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
#include "timer.h"
#include "utilities.h"
#include "variance_reduction.h"
#include "vectors.h"


namespace Impurity
{
	// Function to decide starting t,x,y,z based on input options
	double get_birth_val(const Background::Background& bkg, 
		const int start_opt_int, const double start_val, const double range_min, 
		const double range_max, const BkgFPType bkg_min, 
		const BkgFPType bkg_max)
	{
		// Start at specific point
		double return_start_val {};
		if (start_opt_int == 0)
		{
			 return_start_val = start_val;
		}

		// Start between a user-specified range 
		else if (start_opt_int == 1)
		{
			return_start_val = Random::get(range_min, range_max);
		}

		// Start between the full range of the simulation volume 
		else if (start_opt_int == 2)
		{
			return_start_val = Random::get(static_cast<double>(bkg_min), 
				static_cast<double>(bkg_max));
		}
		
		return return_start_val;
	}

	// Starting t
	double get_birth_t(const Background::Background& bkg, 
		const Options::Options& opts)
	{
		return get_birth_val(bkg, opts.imp_tstart_opt_int(), 
			opts.imp_tstart_val(), opts.imp_trange_min(), opts.imp_trange_max(),
			bkg.get_t_min(), bkg.get_t_max());

		// The random get uses doubles, so need to make sure we are passing
		// double if float is used for BkgFPType
		//return Random::get(static_cast<double>(bkg.get_t_min()), 
		//	static_cast<double>(bkg.get_t_max()));
	}

	// Starting x
	double get_birth_x(const Background::Background& bkg, 
		const Options::Options& opts)
	{
		return get_birth_val(bkg, opts.imp_xstart_opt_int(), 
			opts.imp_xstart_val(), opts.imp_xrange_min(), opts.imp_xrange_max(),
			bkg.get_x_min(), bkg.get_x_max());
	}

	// Starting y
	double get_birth_y(const Background::Background& bkg, 
		const Options::Options& opts)
	{
		return get_birth_val(bkg, opts.imp_ystart_opt_int(), 
			opts.imp_ystart_val(), opts.imp_yrange_min(), opts.imp_yrange_max(),
			bkg.get_y_min(), bkg.get_y_max());
	}

	// Starting z
	double get_birth_z(const Background::Background& bkg, 
		const Options::Options& opts)
	{
		return get_birth_val(bkg, opts.imp_zstart_opt_int(), 
			opts.imp_zstart_val(), opts.imp_zrange_min(), opts.imp_zrange_max(),
			bkg.get_z_min(), bkg.get_z_max());
	}

	std::array<double, 3> get_birth_vXYZ(const Background::Background& bkg, 
		const Options::Options& opts, const double t_imp, const double x_imp, 
		const double y_imp, const double z_imp)
	{
		// If we're running a test we will want to start at predetermined
		// velocities to make sure things work correctly. Ideally there
		// wouldn't be an if statement to favor vectorization, but seeing as
		// this is only happens once per particle it's not a huge deal.
		if (opts.bkg_source_int() == 0)
		{
			// Give particles an initial Y velocity to kick off gyration
			if (opts.test_opt_int() == 0 || opts.test_opt_int() == 1
				|| opts.test_opt_int() == 2 || opts.test_opt_int() == 3)
			{
				return {0.0, 10000.0, 0.0};
			}

			// Curvature drift test requires an initial velocity parallel to
			// the field line, we set that to be 5000 m/s. In this geometry
			// x = R, y = Z and z = phi
			else if (opts.test_opt_int() == 4)
			{
				constexpr double v_par = 30000;
				double v_X {- v_par * std::sin(z_imp)};
				double v_Y {v_par * std::cos(z_imp)};
				return {v_X, v_Y, 0.0};
			}

			// Friction force test case we start at rest so it can accelerate 
			// up to the background velocity
			else if (opts.test_opt_int() == 5)
			{
				return {0.0, 0.0, 0.0};
			}
		}

		// Start at input temperature value
		double start_temp {};
		if (opts.imp_temp_start_opt_int() == 0)
		{
			start_temp = opts.imp_temp_start_val();
		}

		// Start at main ion temperature
		else if (opts.imp_temp_start_opt_int() == 1)
		{
			// Need to get the starting cell indices so we can index the main
			// ion temperature array.
			int tidx {get_nearest_index(bkg.get_times(), 
				static_cast<BkgFPType>(t_imp))};
			int xidx = get_nearest_cell_index(bkg.get_grid_x(),     
				static_cast<BkgFPType>(x_imp));
			int yidx = get_nearest_cell_index(bkg.get_grid_y(),     
				static_cast<BkgFPType>(y_imp));
			int zidx = get_nearest_cell_index(bkg.get_grid_z(),     
				static_cast<BkgFPType>(z_imp));
			
			start_temp = bkg.get_ti()(tidx, xidx, yidx, zidx);
		}

		// Sample from a Maxwellian with mu = sqrt(kT/m) and mean = 0 for an
		// isotropic velocity distribution. Unfortunately throwing away a
		// random number here. T [eV], m [kg]
		const double mu {std::sqrt(start_temp * Constants::ev_to_j 
			/ (opts.imp_mass_amu() * Constants::amu_to_kg))};  // m/s

		double vX {};
		double vY {};
		double vZ {};
		std::tie(vX, vY) = Random::get_two_norm(0.0, mu);
		std::tie(vZ, std::ignore) = Random::get_two_norm(0.0, mu);

		return {vX, vY, vZ};
	}


	int get_birth_charge(const Options::Options& opts)
	{
		return opts.imp_init_charge();
	}

	Impurity create_primary_imp(const Background::Background& bkg, 
		const Options::Options& opts)
	{
		// Starting t,x,y,z for the impurity
		double t_imp = get_birth_t(bkg, opts);
		double x_imp = get_birth_x(bkg, opts);
		double y_imp = get_birth_y(bkg, opts);
		double z_imp = get_birth_z(bkg, opts);

		// Planning to remove this since it's not needed and also tricky
		//auto [X_imp, Y_imp, Z_imp] = opts.mapc2p()(x_imp, y_imp, z_imp);
		double X_imp {0.0};
		double Y_imp {0.0};
		double Z_imp {0.0};

		// Get starting X, Y, Z velocity
		auto [vX_imp, vY_imp, vZ_imp] = get_birth_vXYZ(bkg, opts, t_imp, x_imp, 
			y_imp, z_imp);
		//std::cout << "vX, vY, vZ = " << vX_imp << "\t" << vY_imp << "\t" 
		//	<< vZ_imp << '\n';

		/*
		// Test: Set to true if using a test background plasma
		bool init_test {opts.bkg_source_int() == 0};

		// Test: Slab geometry (x=X, y=Y, z=Z). Testing for gyration (0), 
		// ExB (1) grad-B (2) or polarization (3) drifts, so give particle 
		// initial Y (y) velocity of 2500 m/s to kick off the gyration.
		bool init_slab {(unsigned) < 4u};
		vY_imp = 2500.0 * init_slab * init_test;

		// Test: Cylindrical geometry (x=R, y=Z, z=phi). Testing for curvature 
		// drift (4), so give particle initial velocity of 5000 m/s along the 
		// toroidal field line (the z direction) and 2500 m/s in the Z direction
		// so it gyrates.
		bool init_cyl {(unsigned)x == 4u};
		double vX_imp {-std::sin(z_imp) * 5000.0 * init_cyl * init_test};
		double vY_imp {std::cos(z_imp) * 5000.0 * init_cyl * init_test};
		double vZ_imp {2500.0 * init_cyl * init_test};  // For gyrating
		*/

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

	/*
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
	*/

	// Interpolate the reciprocal basis functions at the impurity location
	std::array<double, 9> interp_recp(const Impurity& imp, 
		const Background::Background& bkg, const int xidx, const int yidx, 
		const int zidx)
	{
		// Get nearest neighbor indices for each direction. These tell us
		// which direction we should interpolate towards, i.e., which
		// rectangle made by the neighboring cell centers our particle
		// is bounded by.
		const int xidx_neighbor {Utilities::get_neighbor_index(imp.get_x(), 
			bkg.get_x(), xidx)};
		const int yidx_neighbor {Utilities::get_neighbor_index(imp.get_y(), 
			bkg.get_y(), yidx)};
		const int zidx_neighbor {Utilities::get_neighbor_index(imp.get_z(), 
			bkg.get_z(), zidx)};

		// x, y, z coordinates of two bounding vertices to interpolate 
		// between. Note these are not grid vertices, but rather are formed
		// by cell center coordinates since that's where B/E are assumed
		// to be defined. It is essentially a cell shifted by dx/2, dy/2
		// and dz/2 if that helps.
		const double x0 {bkg.get_x()[xidx]};
		const double x1 {bkg.get_x()[xidx_neighbor]};
		const double y0 {bkg.get_y()[yidx]};
		const double y1 {bkg.get_y()[yidx_neighbor]};
		const double z0 {bkg.get_z()[zidx]};
		const double z1 {bkg.get_z()[zidx_neighbor]};

		// Array of references to each reciprocal basis vector
		const std::array<std::reference_wrapper<
			const Vectors::Vector3D<BkgFPType>>, 9> recp_basis {
			bkg.get_dxdX(), bkg.get_dxdY(), bkg.get_dxdZ(),
			bkg.get_dydX(), bkg.get_dydY(), bkg.get_dydZ(),
			bkg.get_dzdX(), bkg.get_dzdY(), bkg.get_dzdZ()};

		// Loop through and interpolate each basis vector
		std::array<double, 9> interp_vals {};
		for (int i {}; i < 9; ++i)
		{
			// Values at each vertex, 8 total because it's a rectangle.
			const double v000 {recp_basis[i](xidx, yidx, zidx)};
			const double v100 {recp_basis[i](xidx_neighbor, yidx, zidx)};
			const double v010 {recp_basis[i](xidx, yidx_neighbor, zidx)};
			const double v110 {recp_basis[i](xidx_neighbor, yidx_neighbor, 
				zidx)};
			const double v001 {recp_basis[i](xidx, yidx, zidx_neighbor)};
			const double v101 {recp_basis[i](xidx_neighbor, yidx, 
				zidx_neighbor)};
			const double v011 {recp_basis[i](xidx, yidx_neighbor, 
				zidx_neighbor)};
			const double v111 {recp_basis[i](xidx_neighbor, yidx_neighbor,	
				zidx_neighbor)};

			// Perform interpolation, storing value in our local array
			interp_vals[i] = Utilities::trilinear_interpolate(x0, y0, z0, 
				x1, y1, z1, v000, v100, v010, v110, v001, v101, v011, v111, 
				imp.get_x(), imp.get_y(), imp.get_z());
		}
		return interp_vals;
	}

	bool step(Impurity& imp, const double imp_time_step, 
		const Background::Background& bkg, const Options::Options& opts, 
		const int tidx, int& xidx, int& yidx, int& zidx)
	{
		// If we've run past the time range covered by the background plasma
		// then we're done
		imp.set_t(imp.get_t() + imp_time_step);
		if (imp.get_t() > bkg.get_t_max())
		{
			return false;
		}

		// Step in physical space
        //double dX {imp.get_X() - imp.get_prevX()};
        //double dY {imp.get_Y() - imp.get_prevY()};
        //double dZ {imp.get_Z() - imp.get_prevZ()};

		// Use Boris algorithm to update particle velocity, defined at 
		// half-steps, v_(t-dt/2) --> v_(t+dt/2).
		Boris::update_velocity(imp, bkg, opts, imp_time_step, tidx, xidx, yidx, 
			zidx);

		// Update particle position in physical space after Boris update. This
		// code is actually bad in that it allows a lot of numerical diffusion
		// to enter the update. 
		//imp.set_X(imp.get_X() + imp.get_vX() * imp_time_step);
		//imp.set_Y(imp.get_Y() + imp.get_vY() * imp_time_step);
		//imp.set_Z(imp.get_Z() + imp.get_vZ() * imp_time_step);
        //double dX {imp.get_X() - imp.get_prevX()};
        //double dY {imp.get_Y() - imp.get_prevY()};
        //double dZ {imp.get_Z() - imp.get_prevZ()};
		//double dX {imp.get_vX() * imp_time_step};
		//double dY {imp.get_vY() * imp_time_step};
		//double dZ {imp.get_vZ() * imp_time_step};

		// Calculate the step in computational space using the reciprocal basis
		// vectors. This is Eq. 2.3.6 in Dhaeseleer. Can either use the discrete
		// values for the grid, or interpolate them for a continuous basis
		// vector. Both have their drawbacks though. Discrete is technically
		// only correct at the cell center, and interpolated does not
		// guarantee orthogonality (likely violates it). Pick your poison.
		
		//double dx {bkg.get_dxdX()(xidx, yidx, zidx) * dX
		//	+ bkg.get_dxdY()(xidx, yidx, zidx) * dY
		//	+ bkg.get_dxdZ()(xidx, yidx, zidx) * dZ};
		//double dy {bkg.get_dydX()(xidx, yidx, zidx) * dX
		//	+ bkg.get_dydY()(xidx, yidx, zidx) * dY
		//	+ bkg.get_dydZ()(xidx, yidx, zidx) * dZ};
		//double dz {bkg.get_dzdX()(xidx, yidx, zidx) * dX
		//	+ bkg.get_dzdY()(xidx, yidx, zidx) * dY
		//	+ bkg.get_dzdZ()(xidx, yidx, zidx) * dZ};

		// Interpolate reciprocal basis vector at impurity location
		//std::array<double, 9> int_rec_bas {interp_recp(imp, bkg, xidx, 
		//	yidx, zidx)};
		//double dx {int_rec_bas[0] * dX + int_rec_bas[1] * dY 
		//	+ int_rec_bas[2] * dZ};
		//double dy {int_rec_bas[3] * dX + int_rec_bas[4] * dY 
		//	+ int_rec_bas[5] * dZ};
		//double dz {int_rec_bas[6] * dX + int_rec_bas[7] * dY 
		//	+ int_rec_bas[8] * dZ};

		// Update position in computational space
        //imp.set_x(imp.get_x() + dx);
        //imp.set_y(imp.get_y() + dy);
        //imp.set_z(imp.get_z() + dz);

		// This approach is redundant since vx,vy,vz are pretty much calculated
		// in the dx,dy,dz calculations above. It's just if we want to have 
		// the curvilinear velocity tracked we can do it that way too.
        imp.set_x(imp.get_x() + imp.get_vx() * imp_time_step);
        imp.set_y(imp.get_y() + imp.get_vy() * imp_time_step);
        imp.set_z(imp.get_z() + imp.get_vz() * imp_time_step);

        // Bound checking (move to separate function). 
		// Absorbing at maximum x. xbound_buffer move the BC off the x bound
		// by that much to help avoid some common issues in the background
		// that can happen there, causing impurities to "stick" to the 
		// boundary instead of a proper BC being applied 
		//std::cout << imp.get_x() << "\t" << imp.get_y() << "\t" 
		//	<< imp.get_z() << '\n';
		if ((imp.get_x() + opts.imp_xbound_buffer()) > bkg.get_grid_x().back()) 
		{
			//std::cout << "Absorbed: Max x\n";
			//std::cout << imp.get_x() << "\t" << imp.get_y() << "\t" 
			//	<< imp.get_z() << '\n';
			return false;
		}

		// Minimum x can also be absorbing (0) or mimic a core boundary (1)
		if (opts.min_xbound_type_int() == 0)
		{
			if ((imp.get_x() - opts.imp_xbound_buffer()) < bkg.get_grid_x()[0]) 
			{
				//std::cout << "Absorbed: Min x\n";
				return false;
			}
		}

		else if (opts.min_xbound_type_int() == 1)
		{
			// At minimum x move the particle to a random y,z cell. This is a
			// rough approximation to entering the core and leaving it 
			// somewhere else.
			if ((imp.get_x() - opts.imp_xbound_buffer()) < bkg.get_grid_x()[0])
			{
				//std::cout << "Core BC applied\n";

				imp.set_x(bkg.get_grid_x()[0] + opts.imp_xbound_buffer());
				imp.set_y(Random::get(static_cast<double>(bkg.get_y_min()), 
					static_cast<double>(bkg.get_y_max())));
				imp.set_z(Random::get(static_cast<double>(bkg.get_z_min()), 
					static_cast<double>(bkg.get_z_max())));
			}
		}

		// Separatrix boundary condition. There are two associated values with
		// this option that define the z locations of the X-points. Between
		// those z coordinates a core BC is applied, and outside of them we
		// use an absorbing BC to mimic particle loss to the PFZ. 
		else if (opts.min_xbound_type_int() == 2)
		{
			if ((imp.get_x() - opts.imp_xbound_buffer()) < bkg.get_grid_x()[0])
			{
				// Check if between z extents of X-point (core BC)
				if (imp.get_z() > opts.sep_x_bc_xp_z1() && 
					imp.get_z() < opts.sep_x_bc_xp_z2())
				{
					// Move particle to a random x, y, where z remains in the
					// core range between the X-point z coordinates.
					imp.set_x(bkg.get_grid_x()[0] + opts.imp_xbound_buffer());
					imp.set_y(Random::get(static_cast<double>(bkg.get_y_min()), 
						static_cast<double>(bkg.get_y_max())));
					imp.set_z(Random::get(static_cast<double>
						(opts.sep_x_bc_xp_z1()), 
						static_cast<double>(opts.sep_x_bc_xp_z2())));
				}

				// Otherwise we crossed into the PFZ so count as absorbed.
				return false;
			}
		}

        // Periodic y boundary
        if (imp.get_y() < bkg.get_grid_y()[0])
        {
            imp.set_y(bkg.get_grid_y().back() + (imp.get_y() 
				- bkg.get_grid_y()[0]));
			//std::cout << "Periodic: Min y\n";
        }
        else if (imp.get_y() > bkg.get_grid_y().back())
        {
            imp.set_y(bkg.get_grid_y()[0] + (imp.get_y() 
				- bkg.get_grid_y().back()));
			//std::cout << "Periodic: Max y\n";
        }

		// Absorbing z boundary in SOL, periodic in core
		if (imp.get_x() > opts.lcfs_x())
		{
			if (imp.get_z() < bkg.get_grid_z()[0] || 
				imp.get_z() > bkg.get_grid_z().back()) 
				{
					//std::cout << "Absorbed: Min/max z\n";
					return false;
				}
		}
		else
		{
			if (imp.get_z() < bkg.get_grid_z()[0])
			{
				imp.set_z(bkg.get_grid_z().back() + (imp.get_z() 
					- bkg.get_grid_z()[0]));
				//std::cout << "Periodic: Min z\n";
			}
			else if (imp.get_z() > bkg.get_grid_z().back())
			{
				imp.set_z(bkg.get_grid_z()[0] + (imp.get_z() 
					- bkg.get_grid_z().back()));
				//std::cout << "Periodic: Max z\n";
			}
		}

		return true;
	}

	void record_stats(Statistics& imp_stats, const Impurity& imp, 
		const Background::Background& bkg, const Options::Options& opts, 
		const int tidx, const int xidx, const int yidx, const int zidx, 
		const double imp_time_step)
	{
		// Add one to counts to this location	
		imp_stats.add_counts(tidx, xidx, yidx, zidx, 1);

		// Add total weight to this location. The total weight is the particle
		// weight times the time step. This is Monte Carlo stuff, best read
		// the literature if you don't recognize it.
		imp_stats.add_weights(tidx, xidx, yidx, zidx, 
			static_cast<BkgFPType>(imp.get_weight() * imp_time_step));

		// Add each velocity component to the running sum for this location
		imp_stats.add_vels(tidx, xidx, yidx, zidx, 
			static_cast<BkgFPType>(imp.get_vX()), 
			static_cast<BkgFPType>(imp.get_vY()), 
			static_cast<BkgFPType>(imp.get_vZ()), 
			static_cast<BkgFPType>(imp.get_vx()), 
			static_cast<BkgFPType>(imp.get_vy()), 
			static_cast<BkgFPType>(imp.get_vz()), 
			bkg);

		// Add charge to the running sum for this location
		imp_stats.add_charge(tidx, xidx, yidx, zidx, 
			static_cast<BkgFPType>(imp.get_charge()));

		// Optionally update particle track. This can take a solid chunk of
		// memory so it should only be used with a single particle. It is 
		// mainly used with the test cases so we can show that single particle
		// will have the right ExB drift or that the collision model works,
		// things like that.
		if (opts.save_track_int() == 1) imp_stats.update_track(imp);
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
		const int yidx, const int zidx, const Options::Options& opts,
		std::vector<Impurity>& imps, const std::vector<int>& var_red_counts, 
		Statistics& imp_stats)
	{
		// Update impurity velocity based on Nanbu collision model. Impurity
		// is modified within function. First call is for ions (the false) and
		// second call is for electrons (the true).
		Collisions::nanbu_coll(imp, bkg, tidx, xidx, yidx, zidx, opts, false, 
			imp_time_step, imp_stats);

		// friction_force test case only considers ion collisions to compare
		// against expected flow
		if (opts.test_opt_int() != 5)
		{
			Collisions::nanbu_coll(imp, bkg, tidx, xidx, yidx, zidx, opts, true, 
				imp_time_step, imp_stats);
		}
	}

	void follow_impurity(Impurity& imp, const Background::Background& bkg, 
		Statistics& imp_stats, const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADAS& oa_recomb, int& ioniz_warnings, 
		int& recomb_warnings, std::vector<Impurity>& imps,
		const std::vector<int> var_red_counts, 
		const bool var_red_on, 
		const Options::Options& opts, Timer::Timer& timer)
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
		Boris::update_velocity(imp, bkg, opts, -imp_time_step / 2.0, tidx, xidx, 
			yidx, zidx);

		// Record starting position in statistics arrays
		record_stats(imp_stats, imp, bkg, opts, tidx, xidx, yidx, zidx, 
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
			//auto [fX, fY, fZ] = lorentz_forces(imp, bkg, tidx, xidx, yidx, 
			//	zidx);

			if (var_red_on)
			{
				// Particle splitting based on ionization/recombination 
				// (var_red_split_int == 1). This may split the particle if
				// it is deemed to be in a high importance (e.g., low count) 
				// region according to its ionization or recombination 
				// probability. The split particle is appended to imps to be 
				// followed later.
				if (opts.var_red_split_int() == 1 
					&& opts.imp_iz_recomb_int() > 0)
				{
					VarianceReduction::split_iz_rec(imp, tidx, xidx, yidx, zidx, 
						opts, var_red_counts, imp_time_step, imp_stats, bkg,
						imps, oa_recomb, oa_ioniz);
				}

				// Russian roullete with particle. If particle is in a low
				// importance region (e.g., high counts), then it is killed
				// off with some probability p = var_red_rusrol_prob. If it
				// is not killed off, then the weight increases by 1/p. 
				if (opts.var_red_rusrol_int() > 0) continue_following =
					VarianceReduction::russian_roulette(imp, opts, tidx, xidx, 
					yidx, zidx, var_red_counts, imp_stats);
				if (!continue_following) break;
			}

			// Check for a collision. Since this happens after the Boris update,
			// the velocity stored in imp is v_n+1/2 (and thus prev_v 
			// is v_n-1/2). This is important because collisions happen at
			// full time steps, not the half time steps Boris defines the 
			// velocity at, so the collision step will operate 
			// on v_n = (v_n-1/2 + v_n+1/2) / 2 and then appropriately update
			// v_n+1/2 so the next time step works as expected.
			if (opts.imp_collisions_int() > 0)
			{
				timer.start_coll_timer();
				collision(imp, bkg, imp_time_step, tidx, xidx, yidx, zidx, 
					opts, imps, var_red_counts, imp_stats);
				timer.end_coll_timer();
			}

			// Debugging
			if (debug)
			{
				double R {std::sqrt(imp.get_X()*imp.get_X() 
					+ imp.get_Y()*imp.get_Y() + imp.get_Z()*imp.get_Z())};
				std::cout << "Before step\n";
				std::cout << "  id, tidx, q, t, x, y, z, dt, fX, fY, fZ, " << 
					"vX, vY, vZ, X, Y, Z, R, tidx, xidx, yidx, zidx: \n" 
					<< " " << imp_id << ", " 
					<< tidx << ", " << imp.get_charge() << ", "<< imp.get_t() 
					<< ", " << imp.get_x() << ", " << imp.get_y() << ", " 
					<< imp.get_z() << ", " << imp_time_step << '\n' << "  " 
					<< imp.get_vX() 
					<< ", " << imp.get_vY() << ", " << imp.get_vZ() << ", "
					<< imp.get_X() << ", " << imp.get_Y() << ", " 
					<< imp.get_Z() << ", " << R << '\n' << tidx << " " 
					<< xidx << " " << yidx << " " << zidx << '\n';
			}

			// Perform particle step, update time and tidx. The final location 
			// of the particle could end up out of the grid. That's what the 
			// following code does before recording stats.
			continue_following = step(imp, imp_time_step, bkg, 
				opts, tidx, xidx, yidx, zidx);

			if (debug)
			{
				double R {std::sqrt(imp.get_X()*imp.get_X() 
					+ imp.get_Y()*imp.get_Y() + imp.get_Z()*imp.get_Z())};
				std::cout << "After step\n";
				std::cout << "  id, tidx, q, t, x, y, z, dt, fX, fY, fZ, " << 
					"vX, vY, vZ, X, Y, Z, R, tidx, xidx, yidx, zidx: \n" 
					<< " " << imp_id << ", " 
					<< tidx << ", " << imp.get_charge() << ", "<< imp.get_t() 
					<< ", " << imp.get_x() << ", " << imp.get_y() << ", " 
					<< imp.get_z() << ", " << imp_time_step << '\n' << "  " 
					<< imp.get_vX() 
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
			if (continue_following) record_stats(imp_stats, imp, bkg, opts, 
				tidx, xidx, yidx, zidx, imp_time_step);
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
		Options::Options& opts, Timer::Timer& timer)
	{
		// Variance reduction based on ionization/recombination requires 
		// ionization/recombination to be on (duh)
		if (opts.var_red_split_int() == 1 && opts.imp_iz_recomb_int() == 0)
		{
			std::cout << "Warning! Particle splitting based on "
				<< "ionization/recombination requires imp_ioniz_recomb be set"
				<< "to 'on'. Splitting will be turned off.\n";
			opts.set_var_red_split("off");
		}

		// Variable time step is based on collisions, so those must be on.
		// If not, default to constant time step and alert user.
		if (opts.imp_time_step_opt_int() == 1 
			&& opts.imp_collisions_int() == 0)
		{
			std::cout << "Warning! Impurity time step set to variable, but "
				<< "collisions are off. Changing to constant time step.\n";
			opts.set_imp_time_step_opt("constant");
		}

		// As a safeguard, only save particle tracks if a single particle is
		// used. Otherwise it can easily use way too much memory and blow up
		// the output file.
		if (opts.save_track_int() == 1 && opts.imp_num() > 1)
		{
			std::cout << "Warning! save_track = \"on\" and more than "
				<< "a single particle is tracked (imp_num = " 
				<< opts.imp_num() << "). Only one particle is supported, so "
				<< "imp_num is being decreased to 1.\n";
			opts.set_imp_num(1);
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

		// Likewise declare a reduction operator for combining Timer results
		// of just those processes that exist within parallelized code
		#pragma omp declare reduction(+: Timer::Timer: \
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
			reduction(+: imp_stats, ioniz_warnings, recomb_warnings, timer)
		for (int i = 0; i < opts.imp_num(); ++i)
		{
			// Variance reduction: median mode (var_red_import_int == 0). Based 
			// on seeing if the particle is in a "low count" region, which is 
			// defined as areas with counts less than some fraction of the 
			// median. This section periodically calculates what number of 
			// counts qualifies as a low count region for each thread.
			if ((opts.var_red_split_int() > 0 || opts.var_red_rusrol_int() > 0)
				&& opts.var_red_import_int() == 0 && priv_count > 0)
			{
				// This variable is just an integer saying "every X particles
				// followed, update the variance reduction counts". It needs to
				// be greater than zero.
				int var_red_freq_int {static_cast<int>(
					static_cast<double>(opts.imp_num()) / omp_get_num_threads() 
					* opts.var_red_freq())};
				if (var_red_freq_int == 0) var_red_freq_int = 1;

				// Update variance reduction counts according to the variance
				// reduction frequency (var_red_freq), which is a number 
				// between 0 and 1. So for example, a value of 0.10 would 
				// execute this if statement 10 times during the simulation.
				// The modifier determines what fraction of the median is
				// considered low-count. Lower values = more splitting.
				if (priv_count % var_red_freq_int == 0 && priv_count > 0)
				{
					var_red_counts = 
						VarianceReduction::get_counts(imp_stats, 
						opts.var_red_med_mod());

					// We use this separate boolean (initially false) to
					// control if variance reduction occurs so that the
					// simulation can run through for a bit and learn
					// what a low-count region qualifies as. 
					var_red_control = true;
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
					var_red_counts, var_red_control, opts, timer);
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
		Options::Options& opts, Timer::Timer& timer)
	{
		// Initialize particle statistics vectors, all contained within a
		// Statistics object. Option to control if the three velocity arrays
		// are allocated (to save memory).
		Statistics imp_stats {bkg.get_dim1(), bkg.get_dim2(), 
			bkg.get_dim3(), bkg.get_dim4()};

		// Load OpenADAS data needed for ionization/recombination rates if
		// ionization/recombination is on. Believe it or not, this is the
		// preferred way to do this in C++.
		OpenADAS::OpenADAS oa_ioniz = opts.imp_iz_recomb_int() == 1 
			? OpenADAS::OpenADAS(opts.openadas_root(), opts.openadas_year(), 
				opts.imp_atom_num(), "scd") : OpenADAS::OpenADAS();
		OpenADAS::OpenADAS oa_recomb = opts.imp_iz_recomb_int() == 1 
			? OpenADAS::OpenADAS(opts.openadas_root(), opts.openadas_year(), 
				opts.imp_atom_num(), "acd") : OpenADAS::OpenADAS();

		/*
		if (opts.imp_iz_recomb_int() == 1)
		{
			oa_ioniz = opts.openadas_root(), 
				opts.openadas_year(), opts.imp_atom_num(), "scd";
			oa_recomb = opts.openadas_root(), 
				opts.openadas_year(), opts.imp_atom_num(), "acd";
		}
		*/

		// Execute main particle following loop.
		main_loop(bkg, imp_stats, oa_ioniz, oa_recomb, opts, timer);
	
		// Convert the statistics into meaningful quantities. We are scaling
		// the density by the scale factor.
		std::cout << "Calculating derived quantities...\n";
		std::cout << "  Density...\n";
		imp_stats.calc_density(bkg, opts.imp_num(), 
			opts.imp_source_scale_fact());
		std::cout << "  Velocity...\n";
		imp_stats.calc_vels();
		std::cout << "  Charge...\n";
		imp_stats.calc_charge();
		std::cout << "  Nanbu - s...\n";
		imp_stats.calc_s();
		
		return imp_stats;

	}

}
