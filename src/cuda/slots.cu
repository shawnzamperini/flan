#include <cstdio>

#include "options_device.h"
#include "pcg32.h"
#include "slots.cuh"
#include "slots_device.h"

#include "device_constants.cuh"


namespace Slots
{
	// Function to decide starting t,x,y,z based on input options
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

	
	// Initialize a new particle and return it. Important to pass the rng in
	// as a reference since we change its state each time we pull a random
	// number
	__device__ 
	ParticleInitDevice make_new_particle_cuda(pcg32& rng, 
		const Options::OptionsDevice* opts_d)
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
		p.vx = 0.0;
		p.vy = 0.0;
		p.vz = 0.0;
		p.vX = 0.0;
		p.vY = 5000.0;
		p.vZ = 0.0;
		p.weight = 1.0;
		p.q = opts_d->init_charge;

		return p;
	}

	
	// Kernel to compact dead cells at warp level (what?). 
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
	* @brief GPU kernel to fill slots, replacing dead particles with alive ones.
	*
	* We choose to copy slots_d to the kernel each time this is called since
	* the copy time is << kernel launch time, and copying results in clearer
	* code.
	*/
	__global__ void fill_slots_kernel(SlotsDevice slots_d, 
		int rem_parts, pcg32* rngs_d, 
		const Options::OptionsDevice* opts_d)
	{
		// Global index
		//int i = blockIdx.x * blockDim.x + threadIdx.x;
		//if (i >= slots_d.N) return;
		int idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx >= *slots_d.num_dead) return;
		if (idx >= rem_parts) return;

		// Pick slot index only from those that have been found to be dead.
		int i = slots_d.dead_indices[idx];

		// Skip alive slots and add to count of alive slots (used in progress
		// print out). This atomicAdd isn't a big deal for warp divergence
		// since this thread is already exiting early compared to those that
		// continue past this loop.
		//if (slots_d.state[i] == 0) 
		//{
		//	atomicAdd(slots_d.alive, 1);
		//	return;
		//}

		// Atomically claim a particle index. This prevents the race condition
		// in which multiple threads grab the same index.
		//int idx2 = atomicAdd(slots_d.counter, 1);

		// If we ran out of particles, stop
		//if (idx2 >= rem_parts) return;

		// Create a new particle (device-side initializer), passing in the
		// PCG32 RNG for this thread.
		ParticleInitDevice p = make_new_particle_cuda(rngs_d[i], opts_d);

		// Write particle data
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

		// These get assigned in main_loop
		slots_d.tidx[i] = 0;
		slots_d.xidx[i] = 0;
		slots_d.yidx[i] = 0;
		slots_d.zidx[i] = 0;

		// Mark slot alive
		slots_d.state[i] = 0;

		// Newly alive, count it
		//atomicAdd(slots_d.alive, 1);
	}

	__global__ void all_dead_kernel(SlotsDevice slots_d)
	{
		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= slots_d.N) return;

		// If a single thread has an alive particle (state = 0), then all_dead
		// must be false. This is technically a race condition, but we allow
		// it since only one thread needs to write false for 
		if (slots_d.state[i] == 0) *slots_d.all_dead = false;
	}

	// Replace dead particles with alive ones, as long as remaining particles
	// are greater than zero.
	void fill_slots_gpu(SlotsDevice& slots_d, int& rem_parts, int& alive_slots,
		pcg32* rngs_d, Options::OptionsDevice* opts_d)
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

		// Set counter and number of alive slots to zero
		cudaMemset(slots_d.counter, 0, sizeof(*slots_d.counter));
		cudaMemset(slots_d.alive, 0, sizeof(*slots_d.alive));

		int rem_initial {rem_parts};

		// Call GPU kernel to fill slots.
		int fillBlock = 256;
		int fillGrid  = (num_dead + fillBlock - 1) / fillBlock;
		fill_slots_kernel<<<fillGrid, fillBlock>>>(slots_d, rem_parts, rngs_d, 
			opts_d);

		// Retrieve how many particles were actually filled
		int filled {std::min(rem_initial, num_dead)};

		// Retrieve number of alive slots
		alive_slots = slots_d.N - num_dead + filled;

		// Update remaining particles
		rem_parts = rem_initial - filled;
	}


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

		// Copy all_dead back to host and return
		bool all_dead;
		cudaMemcpy(&all_dead, slots_d.all_dead, sizeof(bool), 
			cudaMemcpyDeviceToHost);

		return all_dead;
	}
}
