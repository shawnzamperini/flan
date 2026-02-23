/**
* @file boris.cpp
*
* @brief Implementation of the Boris algorithm for particle stepping in an
* electromagentic field. See useful introduction and sample code that this
* was inspired by here: https://www.particleincell.com/2011/vxb-rotation/
*/
#include <array>
#include <cmath>

#include "background.h"
#include "constants.h"
#include "impurity.h"
#include "utilities.h"
#include "vectors.h"

namespace Boris
{

	// Calculate cross product and return array
	std::array<double, 3> cross_product(const std::array<double, 3>& a, 
		const std::array<double, 3>& b)
	{
		return 
		{ 
			a[1] * b[2] - a[2] * b[1],
			a[2] * b[0] - a[0] * b[2],
			a[0] * b[1] - a[1] * b[0] 
		};
	}

	// Find nearest neighbor index by essentially seeing which half of the
	// cell a particle is in and returning the index of the neighboring cell
	// closest to that half. Designed to be run once per x,y,z.
	int get_neighbor_index(const double val, 
		const std::vector<BkgFPType>& cell_centers, const int idx)
	{
		// dx > 0 --> (1*2 - 1) = +1
		// dx < 0 --> (0*2 - 1) = -1
		double dx_from_center {val - cell_centers[idx]};
		int side {2 * (dx_from_center > 0.0) - 1};

		// Need to check we aren't at grid edges
		int is_left_edge {(idx == 0)};
		int is_right_edge {(idx == std::ssize(cell_centers) - 1)};

		// This will correctly do -1 --> 1 if we're at the left edge, and
		// 1 --> -1 if we're at the right edge. This effectively means we
		// are using the only neighboring option.
		int offset {
			  side * (1 - is_left_edge - is_right_edge)
			+ 1    * is_left_edge
			+ (-1) * is_right_edge};

		return idx + offset;
	}

	void update_velocity(Impurity::Impurity& imp, 
		const Background::Background& bkg, const double dt, 
		const int tidx, const int xidx, const int yidx, const int zidx)
	{
		// Pull out electric and magnetic field components, and particle
		// velocity components and put into arrays for cleaner code
		/*
		std::array<double, 3> B_local {bkg.get_bX()(tidx, xidx, yidx, zidx), 
			bkg.get_bY()(tidx, xidx, yidx, zidx),
			bkg.get_bZ()(tidx, xidx, yidx, zidx)};
		std::array<double, 3> E_local {bkg.get_eX()(tidx, xidx, yidx, zidx), 
			bkg.get_eY()(tidx, xidx, yidx, zidx),
			bkg.get_eZ()(tidx, xidx, yidx, zidx)};
		*/

		// Get nearest neighbor indices for each direction. These tell us
		// which direction we should interpolate towards, i.e., which
		// rectangle made by the neighboring cell centers our particle
		// is bounded by.
		const int xidx_neighbor {get_neighbor_index(imp.get_x(), bkg.get_x(), 
			xidx)};
		const int yidx_neighbor {get_neighbor_index(imp.get_y(), bkg.get_y(), 
			yidx)};
		const int zidx_neighbor {get_neighbor_index(imp.get_z(), bkg.get_z(), 
			zidx)};

		// Trilinear interpolation to ensure we use a continous B/E in the 
		// algorithm. Without interpolating, the particle will:
		// a) Fail to follow the parallel direction of the field line due to
		//    it being discrete in the code (this is called "numerical
		//    diffusion" and introduces artificial cross-field transport).
		// b) Experience large "kicks" at cell boundaries where B/E is 
		//    discontinuous, which can introduce artificial drifts.
		// 
		// The issue in a) generally will still occur since we are still using
		// an approximation to the field line, but this greatly reduces the
		// numerical diffusion since B is now continous.
		//
		// In the following, we wrap the components into a couple nested arrays
		// so that we can use a single nested loop to iterate through each 
		// B and E component and not have to write the trilinear interpolation
		// loop twice for each vector (or even 6 times, one for each vector
		// component!). It adds a bit of complexity but avoids repeating code.

		// Wrap in convienent array. Need a const reference in reference_wrapper
		// because bkg.get_bX() returns a const ref. Need to match the getter. 
		const std::array<std::reference_wrapper<
			const Vectors::Vector4D<BkgFPType>>, 4> B 
			{bkg.get_bX(), bkg.get_bY(), bkg.get_bZ(), bkg.get_bmag()};
		const std::array<std::reference_wrapper<
			const Vectors::Vector4D<BkgFPType>>, 4> E 
			{bkg.get_eX(), bkg.get_eY(), bkg.get_eZ(), bkg.get_emag()};

		// And put again into single array so we can make an outer loop out
		// of these.
		const std::array outer {&B, &E};

		// Similarly, arrays to hold the output local values
		std::array<double, 4> B_local {};
		std::array<double, 4> E_local {};
		std::array<double, 4>* locals [] = {&B_local, &E_local};

		// Loop through each B/E component, performing a trilinear interpolation
		// to get the local value of each. First loop is to choose B or E.
		for (int j = 0; j < 2; ++j) 
		{ 
			// Index out (B or E), dereference and assign to comps. This is 
			// an object the same as B and E above.
			const auto& comps = *outer[j]; 

			// Get reference to output local values. This is a reference to
			// an array of 3 doubles, B_local or E_local.
			auto& out = *locals[j]; 

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
			//std::cout << "x: " << xidx << "\t" << xidx_neighbor << "\n";
			//std::cout << "y: " << yidx << "\t" << yidx_neighbor << "\n";
			//std::cout << "z: " << zidx << "\t" << zidx_neighbor << "\n";

			// Loop through each component to get interpolated value for each
			for (int i = 0; i < 4; ++i) 
			{
				// Values at each vertex, 8 total because it's a rectangle.
				const double v000 {comps[i](tidx, xidx, yidx, zidx)};
				const double v100 {comps[i](tidx, xidx_neighbor, yidx, zidx)};
				const double v010 {comps[i](tidx, xidx, yidx_neighbor, zidx)};
				const double v110 {comps[i](tidx, xidx_neighbor, yidx_neighbor, 
					zidx)};
				const double v001 {comps[i](tidx, xidx, yidx, zidx_neighbor)};
				const double v101 {comps[i](tidx, xidx_neighbor, yidx, 
					zidx_neighbor)};
				const double v011 {comps[i](tidx, xidx, yidx_neighbor, 
					zidx_neighbor)};
				const double v111 {comps[i](tidx, xidx_neighbor, yidx_neighbor,	
					zidx_neighbor)};

				// Perform interpolation, storing value in our local array
				out[i] = Utilities::trilinear_interpolate(x0, y0, z0, 
					x1, y1, z1, v000, v100, v010, v110, v001, v101, v011, v111, 
					imp.get_x(), imp.get_y(), imp.get_z());
			}

			// We interpolated the magnitude as well. This is because there is
			// no guarantee that the interpolated vector will have the same
			// magnitude (in fact it almost certainly shrink). We are faced
			// with a choice:
			//   a) Leave as-is. B/E are continuous at cell boundaries, but
			//      the magntiude will be off. The size of the gyro-orbit will
			//      be affected as a result.
			//   b) Scale via the interpolated magnitude. B/E will be 
			//      discontinuous at cell boundaries, but definitely not as
			//      discontinous as without interpolation. Gyro-orbit will
			//      be more physically correct.
			// We choose option b), because getting the gyro-orbit correct is
			// central to Flan. Slight discontinuities are not ideal, but the
			// particle spends order of magnitude more time not crossing 
			// a boundary, therefore option b) also minimizes numerical
			// diffusion in that sense as well by ensuring the gyro-orbit is
			// as correct as possible.
			const double out_mag {std::sqrt(out[0] * out[0] + out[1] * out[1] 
				+ out[2] * out[2])};
			//std::cout << "out_mag = " << out_mag << '\n';
			double mask {static_cast<double>(out_mag > 0)};
			for (int i {}; i < 3; ++i)
			{
				// out[3] contains the interpolated magnitude. If out_mag
				// is zero, mask will just multiply by zero to force the vector
				// to zero. It's a divide-by-zero guard.
				out[i] = out[i] * out[3] / std::max(out_mag, 1e-30) * mask;
			}

		}

		/*
		// DEBUG
		for (int j {}; j < 2; ++j)
		{
			if (j == 0) std::cout << "B: ";
			//else std::cout << "E: ";
			for (int i {}; i < 3; ++i)
			{
				if (j == 0) std::cout << B_local[i] << "\t";
				//else std::cout << E_local[i] << "\t";
			}
			//double mag {};
			//if (j == 0) mag = std::sqrt(B_local[0]*B_local[0] 
			//	+ B_local[1]*B_local[1] + B_local[2]*B_local[2]);
			//else mag = std::sqrt(E_local[0]*E_local[0] 
			//	+ E_local[1]*E_local[1] + E_local[2]*E_local[2]);
		}
		*/

		// Normal Boris algorithm continues

		// Store in local variables for conveinence
		double q_m {imp.get_charge() * -Constants::charge_e / imp.get_mass()};
		std::array<double, 3> v {imp.get_vX(), imp.get_vY(), imp.get_vZ()};

		// t vector
		std::array<double, 3> t {};
		for (int i {}; i < 3; ++i) t[i] = q_m * B_local[i] * 0.5 * dt;

		// Magnitude of t, squared
		double tmag2 {t[0]*t[0] + t[1]*t[1] + t[2]*t[2]};

		// s vector
		std::array<double, 3> s {};
		for (int i {}; i < 3; ++i) s[i] = 2 * t[i] / (1.0 + tmag2);

		// v minus
		std::array<double, 3> vminus {};
		for (int i {}; i < 3; ++i) vminus[i] = v[i] + q_m * E_local[i] * 0.5*dt;

		// v prime
		std::array<double, 3> vprime {};
		std::array<double, 3> vminus_cross_t {cross_product(vminus, t)};
		for (int i {}; i < 3; ++i) vprime[i] = vminus[i] + vminus_cross_t[i];

		// v plus
		std::array<double, 3> vplus {};
		std::array<double ,3> vprime_cross_s {cross_product(vprime, s)};
		for (int i {}; i < 3; ++i) vplus[i] = vminus[i] + vprime_cross_s[i];

		// v n+1/2
		// Note we are storing particle velocity at the half time steps before
		// the position. I.e., xi = x(ti), vi = v(ti - dt/2). 
		imp.set_vX(vplus[0] + q_m * E_local[0] * 0.5 * dt);
		imp.set_vY(vplus[1] + q_m * E_local[1] * 0.5 * dt);
		imp.set_vZ(vplus[2] + q_m * E_local[2] * 0.5 * dt);

		/*
		std::cout << tidx << "\t" << xidx << "\t" << yidx << "\t" << zidx << "\t" 
			<< std::scientific << imp.get_t() << "\t" << imp.get_x() 
			<< "\t" << imp.get_y() << "\t" << imp.get_z() << "\t"
			<< imp.get_X() << "\t" << imp.get_Y() << "\t" << imp.get_Z() << "\t"
			<< imp.get_vX() << "\t" << imp.get_vY() << "\t" << imp.get_vZ() << "\t"
			<< B_local[0] << "\t" << B_local[1] << "\t" << B_local[2] <<'\n';
		*/
	}
}


