#include <cstdio>

#include "background_device.h"
#include "slots_device.h"

#include <curand_kernel.h>

#include "boundary.cuh"
#include "device_constants.cuh"

namespace Boundary
{

	// Absorbing boundary condition. Check is value is outside the bound.
	// max_bound = true indicates that bound is a maximum, and that we are
	// checking if value >= bound. Likewise, max_bound = false checks if
	// value =< bound. buffer extends the boundary so that it triggers sooner.
	// This is because sometimes the plasma can behave oddly right at the
	// boundaries, so this can circumvent it.
	__device__ 
	int absorbing_bc(double a, double bound, bool max_bound,
		double buffer = 0.0)
	{
		// Compute both comparisons
		const bool ge {(a + buffer) >= bound};  // for max-bound
		const bool le {(a - buffer) <= bound};  // for min-bound

		// Select the correct one using max_bound as a mask
		const bool dead {max_bound ? ge : le};

		// Returns 0 if alive, 1 if dead
		return static_cast<int>(dead);
	}


	// Periodic boundary condition. This type of BC will never kill a particle,
	// just loop it around to the other side. This function returns the new
	// looped value, or just the same value if a BC was not encountered.
	__device__ 
	double periodic_bc(double a, double amin, double amax)
	{
		const double L = amax - amin;

		const double gt = (a > amax);   // wrap high, subtract L
		const double lt = (a < amin);   // wrap low, add L

		// Return new (or same) position
		return a - gt * L + lt * L;
	}


	// Core boundary condition. The coordinate being checked against is a. b
	// and c are the other two coordinates. If a core condition is met, b is
	// moved to a random value between b_min/b_max, likewise for c. This is
	// essentially like entering the core and then popping out at a random
	// location somewhere else along the core boundary. A residence time
	// could easily be added to this if it was useful.
	__device__
	void core_bc(double& a, double& b, double& c, double a_min,
		double a_max, double buffer, double b_min, double b_max,
		double c_min, double c_max, bool max_bound, curandState& rng_state)
	{
		// Only check the boundary side we were called for
		const bool hit = max_bound ? (a + buffer) > a_max 
			: (a - buffer) < a_min;

		if (hit)
		{
			a = max_bound ? a_max - buffer : a_min + buffer;

			// I leave the commented out CPU code here to point something out.
			// The mt19937 CPU implementation is great, but it doesn't easily
			// port to GPUs. Therefore we decide to rely on the CUDA RNG.
			//b = Random::get(b_min, b_max);
			//c = Random::get(c_min, c_max);
			b = b_min + curand_uniform_double(&rng_state) * (b_max - b_min);
			c = c_min + curand_uniform_double(&rng_state) * (c_max - c_min);
		}
	}


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
		const double lcfs_x)
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
			state = absorbing_bc(slots_d.t[i], Background::d_t_max, true);
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
			new_val = periodic_bc(slots_d.t[i], Background::d_t_min, 
				Background::d_t_max);
			slots_d.t[i] = new_val;
		}

		// --------------------
		// Minimum x boundary
		// --------------------
	
		// Absorbing boundary
		if (min_xbound_type_int == 0)
		{
			state = absorbing_bc(slots_d.x[i], Background::d_x_min, false, 
				imp_xbound_buffer);
			slots_d.state[i] = state;

			// If dead (>0) we can skip the rest of the checks.
			if (state) return;
		}

		// Periodic boundary
		else if (min_xbound_type_int == 1)
		{
			new_val = periodic_bc(slots_d.x[i], Background::d_x_min, 
				Background::d_x_max);
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
			core_bc(xtmp, ytmp, ztmp, Background::d_x_min, Background::d_x_max, 
				imp_xbound_buffer, Background::d_y_min, Background::d_y_max, 
				Background::d_z_min, Background::d_z_max, false, 
				slots_d.rng[i]);

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
			state = absorbing_bc(slots_d.x[i], Background::d_x_max, true, 
				imp_xbound_buffer);
			slots_d.state[i] = state;

			// If dead (>0) we can skip the rest of the checks.
			if (state) return;
		}

		// Periodic boundary
		else if (max_xbound_type_int == 1)
		{
			new_val = periodic_bc(slots_d.x[i], Background::d_x_min, 
				Background::d_x_max);
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
			core_bc(xtmp, ytmp, ztmp, Background::d_x_min, Background::d_x_max, 
				imp_xbound_buffer, Background::d_y_min, Background::d_y_max, 
				Background::d_z_min, Background::d_z_max, true, slots_d.rng[i]);

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
			state = absorbing_bc(slots_d.y[i], Background::d_y_min, false, 
				imp_ybound_buffer);
			slots_d.state[i] = state;

			// If dead (>0) we can skip the rest of the checks.
			if (state) return;
		}

		// Periodic boundary
		else if (min_ybound_type_int == 1)
		{
			new_val = periodic_bc(slots_d.y[i], Background::d_y_min, 
				Background::d_y_max);
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
			core_bc(ytmp, ztmp, xtmp, Background::d_y_min, Background::d_y_max, 
				imp_ybound_buffer, Background::d_z_min, Background::d_z_max, 
				Background::d_x_min, Background::d_x_max, false, 
				slots_d.rng[i]);

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
			state = absorbing_bc(slots_d.y[i], Background::d_y_max, true, 
				imp_ybound_buffer);
			slots_d.state[i] = state;

			// If dead (>0) we can skip the rest of the checks.
			if (state) return;
		}

		// Periodic boundary
		else if (max_ybound_type_int == 1)
		{
			new_val = periodic_bc(slots_d.y[i], Background::d_y_min, 
				Background::d_y_max);
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
			core_bc(ytmp, ztmp, xtmp, Background::d_y_min, Background::d_y_max, 
				imp_ybound_buffer, Background::d_z_min, Background::d_z_max, 
				Background::d_x_min, Background::d_x_max, true, 
				slots_d.rng[i]);

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
			state = absorbing_bc(slots_d.z[i], Background::d_z_min, false, 
				imp_zbound_buffer);
			slots_d.state[i] = state;

			// If dead (>0) we can skip the rest of the checks.
			if (state) return;
		}

		// Periodic boundary
		else if (min_zbound_type_int == 1)
		{
			new_val = periodic_bc(slots_d.z[i], Background::d_z_min, 
				Background::d_z_max);
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
			core_bc(ztmp, xtmp, ytmp, Background::d_z_min, Background::d_z_max, 
				imp_zbound_buffer, Background::d_x_min, Background::d_x_max, 
				Background::d_y_min, Background::d_y_max, false, 
				slots_d.rng[i]);

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
			state = absorbing_bc(slots_d.z[i], Background::d_z_max, true, 
				imp_zbound_buffer);
			slots_d.state[i] = state;

			// If dead (>0) we can skip the rest of the checks.
			if (state) return;
		}

		// Periodic boundary
		else if (max_zbound_type_int == 1)
		{
			new_val = periodic_bc(slots_d.z[i], Background::d_z_min, 
				Background::d_z_max);
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
			core_bc(ztmp, xtmp, ytmp, Background::d_z_min, Background::d_z_max, 
				imp_zbound_buffer, Background::d_x_min, Background::d_x_max, 
				Background::d_y_min, Background::d_y_max, true, 
				slots_d.rng[i]);

			// Write back
			slots_d.x[i] = xtmp;
			slots_d.y[i] = ytmp;
			slots_d.z[i] = ztmp;
		}

	}  //check_bounds_kernel


	void check_bounds_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d, 
		const int tbound_type_int,
		const int min_xbound_type_int, const int max_xbound_type_int,
		const int min_ybound_type_int, const int max_ybound_type_int,
		const int min_zbound_type_int, const int max_zbound_type_int,
		const double imp_xbound_buffer, 
		const double imp_ybound_buffer, 
		const double imp_zbound_buffer, 
		const double lcfs_x)
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
			imp_ybound_buffer, imp_zbound_buffer, lcfs_x);
	}

}
