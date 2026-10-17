/**
* @file  boundary.cu
*
* @brief CUDA routines for handling time and spatial boundary conditions.
*
* GPU version of boundary.cpp. Each time step, every particle is checked
* against the time boundary and the minimum and maximum x, y and z
* boundaries. Each boundary can be absorbing (the particle is killed),
* periodic (the particle is wrapped to the other side) or, for the spatial
* boundaries, core (the particle is put back at the boundary at a random
* location along it). The simulation bounds are read from the __constant__
* variables in device_constants.cuh.
*/
#include <cstdio>

#include "background_device.h"
#include "pcg32.h"
#include "slots_device.h"

#include "boundary.cuh"
#include "device_constants.cuh"

namespace Boundary
{

	/**
	* Checks if a value is outside an absorbing boundary.
	*
	* max_bound = true means bound is a maximum and the check is
	* a + buffer >= bound. max_bound = false means bound is a minimum and
	* the check is a - buffer <= bound. The buffer makes the boundary
	* trigger sooner, since the plasma can behave oddly right at the
	* boundaries. Written without branches so it is SIMD friendly.
	*
	* @param a Value to check (e.g., the particle's x coordinate).
	* @param bound Boundary value.
	* @param max_bound If true, bound is a maximum. Otherwise it is a
	*        minimum.
	* @param buffer Distance the boundary is moved inward. Defaults to 0.
	* @return 1 if the particle is outside the boundary (dead), else 0.
	*/
	__device__
	int absorbing_bc_gpu(double a, double bound, bool max_bound,
		double buffer)
	{
		// Compute both comparisons
		const bool ge {(a + buffer) >= bound};  // for max-bound
		const bool le {(a - buffer) <= bound};  // for min-bound

		// Select the correct one using max_bound as a mask
		const bool dead {max_bound ? ge : le};

		// Returns 0 if alive, 1 if dead
		return static_cast<int>(dead);
	}


	/**
	* Applies a periodic boundary condition.
	*
	* A periodic boundary never kills a particle. A value above amax has
	* the domain length (amax - amin) subtracted, and a value below amin has
	* it added. Only a single wrap is done, so a particle is assumed to
	* move less than one domain length per time step. Written without
	* branches so it is SIMD friendly.
	*
	* @param a Value to check (e.g., the particle's y coordinate).
	* @param amin Minimum of the domain.
	* @param amax Maximum of the domain.
	* @return The wrapped value, or a unchanged if it is inside the domain.
	*/
	__device__
	double periodic_bc_gpu(double a, double amin, double amax)
	{
		const double L = amax - amin;

		const double gt = (a > amax);   // wrap high, subtract L
		const double lt = (a < amin);   // wrap low, add L

		// Return new (or same) position
		return a - gt * L + lt * L;
	}


	/**
	* Applies a core boundary condition.
	*
	* The coordinate being checked is a, and b and c are the other two
	* coordinates. If a is past the boundary (moved inward by a_buffer), a
	* is put back on the buffered boundary, and b and c are moved to
	* uniformly random values between their buffered minimums and
	* maximums. This is like the particle entering the core and popping
	* back out somewhere else along the core boundary. A residence time
	* could be added if it was useful.
	*
	* @param[in,out] a Coordinate checked against the boundary.
	* @param[in,out] b First of the other two coordinates.
	* @param[in,out] c Second of the other two coordinates.
	* @param a_min Minimum of a.
	* @param a_max Maximum of a.
	* @param a_buffer Distance the a boundary is moved inward.
	* @param b_min Minimum of b.
	* @param b_max Maximum of b.
	* @param b_buffer Buffer on each side of the b range.
	* @param c_min Minimum of c.
	* @param c_max Maximum of c.
	* @param c_buffer Buffer on each side of the c range.
	* @param max_bound If true, check the a_max side. Otherwise check the
	*        a_min side.
	* @param rng Random number generator for this particle.
	*/
	__device__
	void core_bc_gpu(double& a, double& b, double& c, double a_min,
		double a_max, double a_buffer, double b_min, double b_max,
		double b_buffer, double c_min, double c_max, double c_buffer, 
		bool max_bound, pcg32& rng)
	{
		// Only check the boundary side we were called for, accounting for 
		// a (symmetric) buffer to avoid artificial edge effects.
		const bool hit = max_bound ? (a + a_buffer) >= a_max 
			: (a - a_buffer) <= a_min;

		if (hit)
		{
			a = max_bound ? a_max - a_buffer : a_min + a_buffer;

			// Particle starts at random b/c locations (accounting for buffers)
			const double b_lo = b_min + b_buffer;
			const double b_hi = b_max - b_buffer;
			b = b_lo + rng.next_double() * (b_hi - b_lo);

			const double c_lo = c_min + c_buffer;
			const double c_hi = c_max - c_buffer;
			c = c_lo + rng.next_double() * (c_hi - c_lo);
		}
	}


	/**
	* Kernel that checks each particle against the time and spatial
	* boundaries and applies the boundary conditions.
	*
	* One thread handles one slot, and dead particles are skipped. The
	* boundaries are checked in the order t, then min/max x, y and z. Each
	* boundary type is an integer: 0 = absorbing, 1 = periodic and
	* 2 = core (x, y and z only). An absorbing boundary sets the slot's
	* state and returns as soon as the particle is dead. For a core
	* boundary in x, y and z are randomized, in y, z and x are randomized,
	* and in z, x and y are randomized.
	*
	* @param slots_d Particle slots on the device. The state, t and
	*        positions are updated in place.
	* @param bkg_d Background plasma on the device. Currently unused.
	* @param tbound_type_int Type of the maximum time boundary (0 or 1).
	* @param min_xbound_type_int Type of the minimum x boundary.
	* @param max_xbound_type_int Type of the maximum x boundary.
	* @param min_ybound_type_int Type of the minimum y boundary.
	* @param max_ybound_type_int Type of the maximum y boundary.
	* @param min_zbound_type_int Type of the minimum z boundary.
	* @param max_zbound_type_int Type of the maximum z boundary.
	* @param imp_xbound_buffer Buffer on the x boundaries.
	* @param imp_ybound_buffer Buffer on the y boundaries.
	* @param imp_zbound_buffer Buffer on the z boundaries.
	* @param lcfs_x x coordinate of the last closed flux surface. Currently
	*        unused.
	* @param rngs_d Random number generators on the device, one per slot.
	*/
	__global__
	void check_bounds_kernel(Slots::SlotsDevice slots_d, 
		const Background::BackgroundDevice bkg_d, 
		const int tbound_type_int,
		const int min_xbound_type_int, const int max_xbound_type_int, 
		const int min_ybound_type_int, const int max_ybound_type_int, 
		const int min_zbound_type_int, const int max_zbound_type_int, 
		const double imp_xbound_buffer, 
		const double imp_ybound_buffer, 
		const double imp_zbound_buffer, 
		const double lcfs_x, pcg32* rngs_d)
	{
		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;

		// Don't try and access beyond the number of slots (segfault)
		if (i >= slots_d.N) return;

		// Skip dead particles
		if (slots_d.state[i] > 0) return;

		// Variables used below
		int state {};
		double new_val {};

		// --------------------
		// Time boundary
		// --------------------

		// Maximum t: Absorbing boundary
		if (tbound_type_int == 0)
		{
			// d_t_max defined in device_constants.cu
			state = absorbing_bc_gpu(slots_d.t[i], d_t_max, true);
			slots_d.state[i] = state;

			// If dead (>0) we can skip the rest of the checks.
			if (state) return;
		}

		// Maximum t: Periodic boundary
		else if (tbound_type_int == 1)
		{
			// We check and overwrite each time here, even though hitting
			// a BC is relatively rare. Despite this, this is still a very
			// fast process since slots.t()[i] is loaded into L1, and 
			// periodic_bc doesn't touch memory loads so the L1 cache never
			// get evicted, so it's a very fast write.
			new_val = periodic_bc_gpu(slots_d.t[i], d_t_min, 
				d_t_max);
			slots_d.t[i] = new_val;
		}

		// --------------------
		// Minimum x boundary
		// --------------------
	
		// Absorbing boundary
		if (min_xbound_type_int == 0)
		{
			state = absorbing_bc_gpu(slots_d.x[i], d_x_min, false, 
				imp_xbound_buffer);
			slots_d.state[i] = state;

			// If dead (>0) we can skip the rest of the checks.
			if (state) return;
		}

		// Periodic boundary
		else if (min_xbound_type_int == 1)
		{
			new_val = periodic_bc_gpu(slots_d.x[i], d_x_min, 
				d_x_max);
			slots_d.x[i] = new_val;
		}

		// Core boundary
		else if (min_xbound_type_int == 2)
		{

			// Temporaries that may get updated in core_bc, which we will
			// write into slots after
			double xtmp {slots_d.x[i]};
			double ytmp {slots_d.y[i]};
			double ztmp {slots_d.z[i]};
			core_bc_gpu(xtmp, ytmp, ztmp, d_x_min, 
				d_x_max, imp_xbound_buffer, d_y_min, 
				d_y_max, imp_ybound_buffer, d_z_min, d_z_max, imp_zbound_buffer,
				false, rngs_d[i]);

			// Write back
			slots_d.x[i] = xtmp;
			slots_d.y[i] = ytmp;
			slots_d.z[i] = ztmp;
		}

		// --------------------
		// Maximum x boundary
		// --------------------
	
		// Absorbing boundary
		if (max_xbound_type_int == 0)
		{
			state = absorbing_bc_gpu(slots_d.x[i], d_x_max, true, 
				imp_xbound_buffer);
			slots_d.state[i] = state;

			// If dead (>0) we can skip the rest of the checks.
			if (state) return;
		}

		// Periodic boundary
		else if (max_xbound_type_int == 1)
		{
			new_val = periodic_bc_gpu(slots_d.x[i], d_x_min, 
				d_x_max);
			slots_d.x[i] = new_val;
		}

		// Core boundary
		else if (max_xbound_type_int == 2)
		{

			// Temporaries that may get updated in core_bc, which we will
			// write into slots after
			double xtmp {slots_d.x[i]};
			double ytmp {slots_d.y[i]};
			double ztmp {slots_d.z[i]};
			core_bc_gpu(xtmp, ytmp, ztmp, d_x_min, 
				d_x_max, imp_xbound_buffer, d_y_min, 
				d_y_max, imp_ybound_buffer,  d_z_min, d_z_max, imp_zbound_buffer,
				true, rngs_d[i]);

			// Write back
			slots_d.x[i] = xtmp;
			slots_d.y[i] = ytmp;
			slots_d.z[i] = ztmp;
		}

		// --------------------
		// Minimum y boundary
		// --------------------
	
		// Absorbing boundary
		if (min_ybound_type_int == 0)
		{
			state = absorbing_bc_gpu(slots_d.y[i], d_y_min, false, 
				imp_ybound_buffer);
			slots_d.state[i] = state;

			// If dead (>0) we can skip the rest of the checks.
			if (state) return;
		}

		// Periodic boundary
		else if (min_ybound_type_int == 1)
		{
			new_val = periodic_bc_gpu(slots_d.y[i], d_y_min, 
				d_y_max);
			slots_d.y[i] = new_val;
		}

		// Core boundary
		else if (min_ybound_type_int == 2)
		{

			// Temporaries that may get updated in core_bc, which we will
			// write into slots after
			double xtmp {slots_d.x[i]};
			double ytmp {slots_d.y[i]};
			double ztmp {slots_d.z[i]};
			core_bc_gpu(ytmp, ztmp, xtmp, d_y_min, 
				d_y_max, imp_ybound_buffer, d_z_min, 
				d_z_max, imp_zbound_buffer, d_x_min, d_x_max, imp_xbound_buffer,
				false, rngs_d[i]);

			// Write back
			slots_d.x[i] = xtmp;
			slots_d.y[i] = ytmp;
			slots_d.z[i] = ztmp;
		}

		// --------------------
		// Maximum y boundary
		// --------------------
	
		// Absorbing boundary
		if (max_ybound_type_int == 0)
		{
			state = absorbing_bc_gpu(slots_d.y[i], d_y_max, true, 
				imp_ybound_buffer);
			slots_d.state[i] = state;

			// If dead (>0) we can skip the rest of the checks.
			if (state) return;
		}

		// Periodic boundary
		else if (max_ybound_type_int == 1)
		{
			new_val = periodic_bc_gpu(slots_d.y[i], d_y_min, 
				d_y_max);
			slots_d.y[i] = new_val;
		}

		// Core boundary
		else if (max_ybound_type_int == 2)
		{

			// Temporaries that may get updated in core_bc, which we will
			// write into slots after
			double xtmp {slots_d.x[i]};
			double ytmp {slots_d.y[i]};
			double ztmp {slots_d.z[i]};
			core_bc_gpu(ytmp, ztmp, xtmp, d_y_min, 
				d_y_max, imp_ybound_buffer, d_z_min, 
				d_z_max, imp_zbound_buffer, d_x_min, d_x_max, imp_xbound_buffer,
				true, rngs_d[i]);

			// Write back
			slots_d.x[i] = xtmp;
			slots_d.y[i] = ytmp;
			slots_d.z[i] = ztmp;
		}

		// --------------------
		// Minimum z boundary
		// --------------------
	
		// Absorbing boundary
		if (min_zbound_type_int == 0)
		{
			state = absorbing_bc_gpu(slots_d.z[i], d_z_min, false, 
				imp_zbound_buffer);
			slots_d.state[i] = state;

			// If dead (>0) we can skip the rest of the checks.
			if (state) return;
		}

		// Periodic boundary
		else if (min_zbound_type_int == 1)
		{
			new_val = periodic_bc_gpu(slots_d.z[i], d_z_min, 
				d_z_max);
			slots_d.z[i] = new_val;
		}

		// Core boundary
		else if (min_zbound_type_int == 2)
		{

			// Temporaries that may get updated in core_bc, which we will
			// write into slots after
			double xtmp {slots_d.x[i]};
			double ytmp {slots_d.y[i]};
			double ztmp {slots_d.z[i]};
			core_bc_gpu(ztmp, xtmp, ytmp, d_z_min, 
				d_z_max, imp_zbound_buffer, d_x_min, 
				d_x_max, imp_xbound_buffer, d_y_min, d_y_max, imp_ybound_buffer,
				false, rngs_d[i]);

			// Write back
			slots_d.x[i] = xtmp;
			slots_d.y[i] = ytmp;
			slots_d.z[i] = ztmp;
		}

		// --------------------
		// Maximum z boundary
		// --------------------
	
		// Absorbing boundary
		if (max_zbound_type_int == 0)
		{
			state = absorbing_bc_gpu(slots_d.z[i], d_z_max, true, 
				imp_zbound_buffer);
			slots_d.state[i] = state;

			// If dead (>0) we can skip the rest of the checks.
			if (state) return;
		}

		// Periodic boundary
		else if (max_zbound_type_int == 1)
		{
			new_val = periodic_bc_gpu(slots_d.z[i], d_z_min, 
				d_z_max);
			slots_d.z[i] = new_val;
		}

		// Core boundary
		else if (max_zbound_type_int == 2)
		{

			// Temporaries that may get updated in core_bc, which we will
			// write into slots after
			double xtmp {slots_d.x[i]};
			double ytmp {slots_d.y[i]};
			double ztmp {slots_d.z[i]};
			core_bc_gpu(ztmp, xtmp, ytmp, d_z_min, 
				d_z_max, imp_zbound_buffer, d_x_min, 
				d_x_max, imp_xbound_buffer, d_y_min, d_y_max, imp_ybound_buffer, 
				true, rngs_d[i]);

			// Write back
			slots_d.x[i] = xtmp;
			slots_d.y[i] = ytmp;
			slots_d.z[i] = ztmp;
		}

	}  //check_bounds_kernel


	/**
	* Launches check_bounds_kernel to apply the time and spatial boundary
	* conditions to every particle on the GPU.
	*
	* Uses blocks of 256 threads, one thread per slot. The launch is
	* asynchronous, except in DEBUG builds, where it waits for the kernel to
	* finish and prints any error.
	*
	* @param slots_d Particle slots on the device.
	* @param bkg_d Background plasma on the device.
	* @param tbound_type_int Type of the maximum time boundary.
	* @param min_xbound_type_int Type of the minimum x boundary.
	* @param max_xbound_type_int Type of the maximum x boundary.
	* @param min_ybound_type_int Type of the minimum y boundary.
	* @param max_ybound_type_int Type of the maximum y boundary.
	* @param min_zbound_type_int Type of the minimum z boundary.
	* @param max_zbound_type_int Type of the maximum z boundary.
	* @param imp_xbound_buffer Buffer on the x boundaries.
	* @param imp_ybound_buffer Buffer on the y boundaries.
	* @param imp_zbound_buffer Buffer on the z boundaries.
	* @param lcfs_x x coordinate of the last closed flux surface.
	* @param rngs_d Random number generators on the device, one per slot.
	*/
	void check_bounds_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d, 
		const int tbound_type_int,
		const int min_xbound_type_int, const int max_xbound_type_int,
		const int min_ybound_type_int, const int max_ybound_type_int,
		const int min_zbound_type_int, const int max_zbound_type_int,
		const double imp_xbound_buffer, 
		const double imp_ybound_buffer, 
		const double imp_zbound_buffer, 
		const double lcfs_x, pcg32* rngs_d)
	{

		// Block and grid size
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;

		// Call kernel to check if particles have encountered boundary condition
		check_bounds_kernel<<<gridSize, blockSize>>>(slots_d, bkg_d, 
			tbound_type_int,
			min_xbound_type_int, max_xbound_type_int, 
			min_ybound_type_int, max_ybound_type_int, 
			min_zbound_type_int, max_zbound_type_int, imp_xbound_buffer, 
			imp_ybound_buffer, imp_zbound_buffer, lcfs_x, rngs_d);

#ifdef DEBUG
		// Check for errors
		cudaError_t err {cudaDeviceSynchronize()};
		if (err != cudaSuccess)
			printf("check_bounds_kernel error: %s\n", cudaGetErrorString(err));
#endif
	}

}
