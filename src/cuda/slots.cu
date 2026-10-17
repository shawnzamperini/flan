/**
* @file  slots.cu
*
* @brief CUDA routines for filling particle slots with new particles and
*        checking if every slot is dead.
*
* GPU version of the slot routines in slots.cpp. Each time step, the
* indices of the dead slots are gathered into a compact list with a
* warp-level stream compaction, and then one thread per dead slot creates a
* new particle there, as long as there are particles left to launch. New
* particles start according to the imp_*start_opt and imp_temp_start_opt
* options.
*/
#include <cstdio>

#include "background_device.h"
#include "options_device.h"
#include "pcg32.h"
#include "slots.cuh"
#include "slots_device.h"

#include "device_constants.cuh"
#include "indices.cuh"
#include "interpolate.cuh"


namespace Slots
{

	/**
	* Simple struct to hold the three components of a vector.
	*/
	struct Vec3d
	{
		double x;
		double y;
		double z;
	};

	/**
	* Picks the starting value of one coordinate (t, x, y or z) of a new
	* particle.
	*
	* start_opt_int is the control integer of the imp_*start_opt option
	* deciding where the particle starts:
	*   0: At start_val.
	*   1: Uniformly random between range_min and range_max.
	*   2: Uniformly random over the full simulation volume, between bkg_min
	*      and bkg_max.
	* Any other value returns 0.
	*
	* @param start_opt_int Control integer of the imp_*start_opt option.
	* @param start_val Starting value used when start_opt_int is 0.
	* @param range_min Lower end of the range used when start_opt_int is 1.
	* @param range_max Upper end of the range used when start_opt_int is 1.
	* @param bkg_min Minimum of the simulation volume, used when start_opt_int
	*        is 2.
	* @param bkg_max Maximum of the simulation volume, used when start_opt_int
	*        is 2.
	* @param rng Random number generator for this thread.
	* @return Starting value of the coordinate.
	*/
	__device__
	double get_birth_val_cuda(const int start_opt_int, const double start_val, 
		const double range_min, const double range_max, const double bkg_min, 
		const double bkg_max, pcg32& rng)
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
			return_start_val = range_min + (range_max - range_min) 
				* rng.next_double();
		}

		// Start between the full range of the simulation volume 
		else if (start_opt_int == 2)
		{
			return_start_val = bkg_min + (bkg_max - bkg_min) 
				* rng.next_double();
		}
		
		return return_start_val;
	}


	/**
	* Picks the starting Cartesian velocity of a new particle.
	*
	* For the test backgrounds (bkg_source = "test") the velocity is fixed
	* by the test case: 10,000 m/s in Y for the gyrate, ExB, grad-B and
	* polarization tests, 30,000 m/s along the field line (in the X-Y plane,
	* at angle z) for the curvature test, and at rest for the friction test.
	*
	* Otherwise each component is sampled from a Maxwellian with zero mean
	* and sigma = sqrt(kT/m), giving an isotropic velocity distribution. The
	* temperature is imp_temp_start_val if imp_temp_start_opt is
	* "single_value", or the main ion temperature interpolated (trilinear
	* in space, linear in time) at the particle location if it is
	* "main_ion". The interpolated ion temperature is held at or above
	* min_ti, since it can go negative near the grid edges.
	*
	* @param bkg_d Background plasma on the device, used for the ion
	*        temperature.
	* @param t Starting time of the particle (s).
	* @param x Starting computational x coordinate.
	* @param y Starting computational y coordinate.
	* @param z Starting computational z coordinate. In the curvature test,
	*        z is the toroidal angle.
	* @param tidx Index of the nearest background time frame.
	* @param xidx x index of the grid cell containing the particle.
	* @param yidx y index of the grid cell containing the particle.
	* @param zidx z index of the grid cell containing the particle.
	* @param rng Random number generator for this thread.
	* @param opts_d Simulation options on the device.
	* @return Starting velocity {vX, vY, vZ} in Cartesian space (m/s).
	*/
	__device__
	Vec3d get_birth_vXYZ(const Background::BackgroundDevice bkg_d, 
		const double t, const double x, const double y, const double z, 
		const int tidx, const int xidx, const int yidx, const int zidx, 
		pcg32& rng, const Options::OptionsDevice* opts_d)
	{
		// If we're running a test we will want to start at predetermined
		// velocities to make sure things work correctly. Ideally there
		// wouldn't be an if statement to favor vectorization, but seeing as
		// this is only happens once per particle it's not a huge deal.
		if (opts_d->bkg_source_int == 0)
		{
			// Give particles an initial Y velocity to kick off gyration
			if (opts_d->test_opt_int == 0 || opts_d->test_opt_int == 1
				|| opts_d->test_opt_int == 2 || opts_d->test_opt_int == 3)
			{
				return {0.0, 10000.0, 0.0};
			}

			// Curvature drift test requires an initial velocity parallel to
			// the field line, we set that to be 30,000 m/s. In this geometry
			// x = R, y = Z and z = phi
			else if (opts_d->test_opt_int == 4)
			{
				constexpr double v_par = 30000;
				double v_X {- v_par * sin(z)};
				double v_Y {v_par * cos(z)};
				return {v_X, v_Y, 0.0};
			}

			// Friction force test case we start at rest so it can accelerate 
			// up to the background velocity
			else if (opts_d->test_opt_int == 5)
			{
				return {0.0, 0.0, 0.0};
			}
		}

		// Start at input temperature value
		double start_temp {};
		if (opts_d->imp_temp_start_opt_int == 0)
		{
			start_temp = opts_d->imp_temp_start_val;
		}

		// Start at main ion temperature, interpolating for it
		// Get interpolation stencil at this location for 4D arrays
		else if (opts_d->imp_temp_start_opt_int == 1)
		{
			Interpolate::InterpolationStencil stencil_4d
				{Interpolate::build_stencil_4d(bkg_d, tidx, xidx, yidx, zidx, 
				t, x, y, z)};

			start_temp = Interpolate::interpolate_field_4d(bkg_d.ti, 
				stencil_4d, t, x, y, z);

			// Limit to minimum value
			start_temp = fmax(start_temp, opts_d->min_ti);
		}

		// Sample from a Maxwellian with sigma = sqrt(kT/m) and mean = 0 for an
		// isotropic velocity distribution. T [eV], m [kg]
		const double sigma {sqrt(start_temp * Constants::ev_to_j 
			/ (opts_d->imp_mass_amu * Constants::amu_to_kg))};  // m/s

		double vX {rng.normal(0.0, sigma)};
		double vY {rng.normal(0.0, sigma)};
		double vZ {rng.normal(0.0, sigma)};

		return {vX, vY, vZ};
	}

	
	/**
	* Creates a new particle based on the input options.
	*
	* The starting t, x, y and z come from get_birth_val_cuda, the indices
	* of the nearest time frame and the containing grid cell are found, and
	* the starting Cartesian velocity comes from get_birth_vXYZ. The
	* computational velocities are left at zero, since they are set on the
	* first Boris update. The weight is 1 / imp_num and the charge is
	* init_charge.
	*
	* @param rng_d Random number generator of the slot the particle goes
	*        into.
	* @param opts_d Simulation options on the device.
	* @param bkg_d Background plasma on the device.
	* @return The new particle.
	*/
	__device__
	ParticleInitDevice make_new_particle_cuda(pcg32& rng, 
		const Options::OptionsDevice* opts_d, 
		const Background::BackgroundDevice& bkg_d)
	{

		ParticleInitDevice p;

		// Create a new particle based on input options.
		// d_x_min/max, etc. are defined in device_constants.cuh
		p.t = get_birth_val_cuda(opts_d->tstart_opt_int, opts_d->tstart_val, 
			opts_d->trange_min, opts_d->trange_max, d_t_min, d_t_max, rng);
		p.x = get_birth_val_cuda(opts_d->xstart_opt_int, opts_d->xstart_val, 
			opts_d->xrange_min, opts_d->xrange_max, d_x_min, d_x_max, rng);
		p.y = get_birth_val_cuda(opts_d->ystart_opt_int, opts_d->ystart_val, 
			opts_d->yrange_min, opts_d->yrange_max, d_y_min, d_y_max, rng);
		p.z = get_birth_val_cuda(opts_d->zstart_opt_int, opts_d->zstart_val, 
			opts_d->zrange_min, opts_d->zrange_max, d_z_min, d_z_max, rng);

		// Get indices, same process as in impurity_transport.cu
		p.tidx = Indices::get_nearest_index_cuda(d_t, bkg_d.tdim, p.t);
		p.xidx = Indices::get_nearest_cell_index_cuda(d_grid_x, bkg_d.xdim+1, 
			p.x);
		p.yidx = Indices::get_nearest_cell_index_cuda(d_grid_y, bkg_d.ydim+1, 
			p.y);
		p.zidx = Indices::get_nearest_cell_index_cuda(d_grid_z, bkg_d.zdim+1, 
			p.z);

		// Get starting vX, vY, vZ into struct
		Vec3d vXYZ {get_birth_vXYZ(bkg_d, p.t, p.x, p.y, p.z, p.tidx, p.xidx, 
			p.yidx, p.zidx, rng, opts_d)};
		p.vX = vXYZ.x;
		p.vY = vXYZ.y;
		p.vZ = vXYZ.z;

		// These get set on the first call to Boris::update_velocity 
		p.vx = 0.0;
		p.vy = 0.0;
		p.vz = 0.0;

		// Particle weights are in units of particles/s, meant to represent a
		// given source of particles at this location. We need to normalize
		// by number of particle launched as well.
		p.weight = 1.0 / opts_d->imp_num;
		p.q = opts_d->init_charge;

		return p;
	}

	
	/**
	* Kernel that gathers the indices of every dead slot into a compact
	* list.
	*
	* This is a warp-level stream compaction. Each warp votes on which of
	* its lanes hold a dead slot (__ballot_sync), lane 0 reserves room for
	* that many indices in dead_indices with a single atomicAdd on
	* num_dead, and each dead lane then writes its slot index at its offset
	* within the warp's reserved block. Using one atomic per warp instead of
	* one per thread keeps contention low. The order of the indices between
	* warps is not deterministic.
	*
	* @param state State of each slot, where 0 is alive and anything else
	*        is dead.
	* @param dead_indices Output list of dead slot indices. Must be able to
	*        hold N indices.
	* @param num_dead Device counter of dead slots. Must be zeroed before
	*        the launch. Holds the number of dead slots afterward.
	* @param N Number of slots.
	*/
	__global__ void compact_dead_slots(const int* __restrict__ state,
		int* __restrict__ dead_indices, int* __restrict__ num_dead, int N)
	{
		// Global thread index
		int tid  = blockIdx.x * blockDim.x + threadIdx.x;

		// Thread position within its warp (lane within warp)
		int lane = threadIdx.x & 31;  // This is essentially thread % 32

		// Active mask for this warp. Tells us which lanes are active.
		// This returns a bitwise mask in a 32-bit integer. For
		// example, if mask = 0b00000000_00000000_00000000_11101111, this 
		// means lanes 0, 1, 2, 3, 5, 6, 7 are active and 4, 8-31 are inactive.
		// Really just relevant for the final warp if slots doesn't equally
		// divide by 32. We need this for the warp-level intrinsics below.
		unsigned int mask {__activemask()};

		// Is this thread in range?
		bool in_range {(tid < N)};

		// Is this slot dead (state > 0)?
		bool is_dead {false};
		if (in_range) is_dead = (state[tid] != 0);

		// Warp-wide vote: which lanes are dead? Similar to activemask, this
		// returns a 32-bit bitwise mask indicating the lanes where
		// is_dead=true. It only evaluates this for active lanes, as set in
		// mask. For example, say our active lanes are:
		//       mask = 0b00000000_00000000_00000000_11101111
		//                                           TFT FFTT  <- is_dead values
		// Then we'd get:
		//  dead_mask = 0b00000000_00000000_00000000_10100011
		unsigned int dead_mask = __ballot_sync(mask, is_dead);

		// Number of dead lanes in this warp (popc = population count). This 
		// just counts the number of 1's in dead_mask. Following the above
		// example, we'd get 4 for this.
		int warp_dead_count = __popc(dead_mask);

		// No dead slots in this warp, nothing to do
		if (warp_dead_count == 0) return;

		// Compute this lane's offset among dead lanes in the warp
		// Mask of dead lanes with index < lane
		// First, 1u << lane creates a 32-bit mask with a single bit at position
		// lane. So if lane = 16, then
		//   1u << lane = 0b00000000_00000001_00000000_00000000
		// Then subtracting 1 from this turns it into a mask with all bits
		// below our bit equal to 1;
		//   (1u << lane - 1) = 0b00000000_00000000_11111111_11111111
		// Then we perform a bitwise & with dead_mask which keeps only the 
		// dead bits below this lane. For our example from above,
		//  dead_mask = 0b00000000_00000000_00000000_10100011
		//              0b00000000_00000000_11111111_11111111
		//  lane_mask = 0b00000000_00000000_00000000_10100011
		// Finally, lane_offset is the sum of all the 1's, lane_offset = 4.
		// This is used in calculating the index in dead_indices that we will 
		// place the dead slot index in, done below.
		unsigned int lane_mask = dead_mask & ((1u << lane) - 1);
		int lane_offset = __popc(lane_mask);

		// Let lane 0 reserve space in the global dead_indices array
		// atomicAdd returns the original value of num_dead, which we assign
		// to warp_base. Continuing our example where warp_dead_count=4, if
		// we just say num_dead was 10 going into this, then after
		// num_dead = 14 and warp_base = 10. Therefore when we get down to
		// dead_indices, we will know that this warp "owns" (will only be
		// indexing into) indices 10, 11, 12, 13, where each lane knows which
		// index it will do via lane_offset.
		int warp_base = 0;
		if (lane == 0) warp_base = atomicAdd(num_dead, warp_dead_count);

		// Broadcast warp_base to all lanes so that every lane in this warp 
		// knows where in dead_indices to start from
		warp_base = __shfl_sync(mask, warp_base, 0);

		// If this lane is dead, write its global index into the compacted array
		if (is_dead) dead_indices[warp_base + lane_offset] = tid;
	}


	/**
	* Kernel that fills dead slots with new particles.
	*
	* One thread handles one entry of the dead slot list from
	* compact_dead_slots, so only the first min(num_dead, rem_parts) dead
	* slots are filled. Each creates a new particle with
	* make_new_particle_cuda and writes it into the slot, marking it alive.
	* The cell indices are set to 0 here and assigned afterward by
	* find_containing_cell_gpu.
	*
	* slots_d is passed by value each time since the copy time is much less
	* than the kernel launch time, and it makes for clearer code.
	*
	* @param slots_d Particle slots on the device. The dead slots are
	*        overwritten in place.
	* @param bkg_d Background plasma on the device.
	* @param rem_parts Number of particles left to launch.
	* @param rngs_d Random number generators on the device.
	* @param opts_d Simulation options on the device.
	*/
	__global__ void fill_slots_kernel(SlotsDevice slots_d, 
		Background::BackgroundDevice bkg_d, int rem_parts, pcg32* rngs_d, 
		const Options::OptionsDevice* opts_d)
	{
		// Global index
		int idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx >= *slots_d.num_dead) return;
		if (idx >= rem_parts) return;

		// Pick slot index only from those that have been found to be dead.
		int i = slots_d.dead_indices[idx];

		// Create a new particle (device-side initializer)
		ParticleInitDevice p = make_new_particle_cuda(rngs_d[i], opts_d, bkg_d);

		// Write particle data
		slots_d.pt[i] = 0.0;
		slots_d.t[i]  = p.t;

		slots_d.x[i]  = p.x;
		slots_d.y[i]  = p.y;
		slots_d.z[i]  = p.z;

		slots_d.vx[i] = p.vx;
		slots_d.vy[i] = p.vy;
		slots_d.vz[i] = p.vz;

		slots_d.vX[i] = p.vX;
		slots_d.vY[i] = p.vY;
		slots_d.vZ[i] = p.vZ;
		
		slots_d.weight[i] = p.weight;

		slots_d.q[i] = p.q;

		// These get assigned in find_containing_cell_gpu
		slots_d.tidx[i] = 0;
		slots_d.xidx[i] = 0;
		slots_d.yidx[i] = 0;
		slots_d.zidx[i] = 0;

		// Mark slot alive
		slots_d.state[i] = 0;
	}

	/**
	* Kernel that sets all_dead to false if any slot holds an alive
	* particle.
	*
	* all_dead must be set to true before the launch. One thread checks
	* one slot.
	*
	* @param slots_d Particle slots on the device.
	*/
	__global__ void all_dead_kernel(SlotsDevice slots_d)
	{
		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= slots_d.N) return;

		// If a single thread has an alive particle (state = 0), then all_dead
		// must be false. This is technically a race condition, but we allow
		// it since every thread that writes is writing the same value (false).
		if (slots_d.state[i] == 0) *slots_d.all_dead = false;
	}

	/**
	* Replaces dead particles with new ones, as long as there are particles
	* left to launch.
	*
	* Runs compact_dead_slots to gather the dead slot indices, copies the
	* number of dead slots back to the host (which waits for the kernel to
	* finish), and then launches fill_slots_kernel with one thread per dead
	* slot. Both use blocks of 256 threads. The fill kernel launch is
	* asynchronous, except in DEBUG builds, where it waits for the kernel to
	* finish and prints any error. The number of filled slots is
	* min(rem_parts, num_dead), so the counts are updated on the host
	* without waiting for it.
	*
	* @param slots_d Particle slots on the device.
	* @param rem_parts Number of particles left to launch. Decreased by the
	*        number of slots filled.
	* @param[out] alive_slots Number of alive slots after filling.
	* @param rngs_d Random number generators on the device.
	* @param opts_d Simulation options on the device.
	* @param bkg_d Background plasma on the device.
	*/
	void fill_slots_gpu(SlotsDevice& slots_d, int& rem_parts, int& alive_slots,
		pcg32* rngs_d, Options::OptionsDevice* opts_d,
		const Background::BackgroundDevice& bkg_d)
	{

		// Each slot is assigned to a specific GPU where its data lies
		cudaSetDevice(slots_d.device_id);

		// Reset number of dead slots to zero
		cudaMemset(slots_d.num_dead, 0, sizeof(*slots_d.num_dead));

		// Fill in dead_indices with just the slots indices that have dead
		// particles.
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;
		compact_dead_slots<<<gridSize, blockSize>>>(slots_d.state,
			slots_d.dead_indices, slots_d.num_dead, slots_d.N);

		// Read back num_dead
		int num_dead = 0;
		cudaMemcpy(&num_dead, slots_d.num_dead, sizeof(*slots_d.num_dead), 
			cudaMemcpyDeviceToHost);

		int rem_initial {rem_parts};

		// Call GPU kernel to fill slots.
		int fillBlock = 256;
		int fillGrid  = (num_dead + fillBlock - 1) / fillBlock;
		fill_slots_kernel<<<fillGrid, fillBlock>>>(slots_d, bkg_d, rem_parts,
			rngs_d, opts_d);

#ifdef DEBUG
		// Check for errors. This also catches errors from compact_dead_slots.
		cudaError_t err {cudaDeviceSynchronize()};
		if (err != cudaSuccess)
			printf("fill_slots_kernel error: %s\n", cudaGetErrorString(err));
#endif

		// Retrieve how many particles were actually filled
		int filled {std::min(rem_initial, num_dead)};

		// Retrieve number of alive slots
		alive_slots = slots_d.N - num_dead + filled;

		// Update remaining particles
		rem_parts = rem_initial - filled;
	}


	/**
	* Checks if every slot holds a dead particle.
	*
	* Sets all_dead to true on the device, launches all_dead_kernel (blocks
	* of 256 threads, one thread per slot) and copies the result back,
	* which waits for the kernel to finish. In DEBUG builds, any kernel
	* error is printed.
	*
	* @param slots_d Particle slots on the device.
	* @return true if every slot is dead, false if not.
	*/
	bool all_dead_gpu(SlotsDevice& slots_d)
	{

		// Each slot is assigned to a specific GPU where its data lies
		cudaSetDevice(slots_d.device_id);

		// Must set all_dead to true each time. Think of it like a question
		// that the device responds to.
		//   Host: Everything is dead? (all_dead = true)
		//   Device: No (all_dead = false)
		//   Device: Yes (all_dead untouched, stays true)
		// As we come back again during the next time step, we have to ask the
		// question again.
		bool init = true;
		cudaMemcpy(slots_d.all_dead, &init, sizeof(bool), 
			cudaMemcpyHostToDevice);

		// Call GPU kernel
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;
		all_dead_kernel<<<gridSize, blockSize>>>(slots_d);

#ifdef DEBUG
		// Check for errors
		cudaError_t err {cudaDeviceSynchronize()};
		if (err != cudaSuccess)
			printf("all_dead_kernel error: %s\n", cudaGetErrorString(err));
#endif

		// Copy all_dead back to host and return
		bool all_dead;
		cudaMemcpy(&all_dead, slots_d.all_dead, sizeof(bool), 
			cudaMemcpyDeviceToHost);

		return all_dead;
	}
}
