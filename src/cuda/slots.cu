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

	/**
	* @brief GPU kernel to fill slots, replacing dead particles with alive ones.
	*
	* We choose to copy slots_d to the kernel each time this is called since
	* the copy time is << kernel launch time, and copying results in clearer
	* code.
	*/
	__global__ void fill_slots_kernel(SlotsDevice slots_d, int* counter, 
		int* alive_counter, int rem_parts, pcg32* rngs_d, 
		const Options::OptionsDevice* opts_d)
	{
		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= slots_d.N) return;

		// Skip alive slots and add to count of alive slots (used in progress
		// print out). This atomicAdd isn't a big deal for warp divergence
		// since this thread is already exiting early compared to those that
		// continue past this loop.
		if (slots_d.state[i] == 0) 
		{
			atomicAdd(alive_counter, 1);
			return;
		}

		// Atomically claim a particle index. This prevents the race condition
		// in which multiple threads grab the same index.
		int idx = atomicAdd(counter, 1);

		// If we ran out of particles, stop
		if (idx >= rem_parts) return;

		// Create a new particle (device-side initializer), passing in the
		// PCG32 RNG for this thread.
		ParticleInitDevice p = make_new_particle_cuda(rngs_d[i], opts_d);

		// Write particle data
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

		// Mark slot alive
		slots_d.state[i] = 0;

		// Newly alive, count it
		atomicAdd(alive_counter, 1);
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
		// No more particles left, don't fill
		if (rem_parts <= 0) return;

		// Each slot is assigned to a specific GPU where its data lies
		cudaSetDevice(slots_d.device_id);

		// Device counter
		int* d_counter;
		cudaMalloc(&d_counter, sizeof(int));

		// Set counter to zero
		int zero = 0;
		cudaMemcpy(d_counter, &zero, sizeof(int), cudaMemcpyHostToDevice);

		// Counter for number of alive slots
		int* d_alive;
		cudaMalloc(&d_alive, sizeof(int));
		cudaMemset(d_alive, 0, sizeof(int));

		// Call GPU kernel to fill slots.
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;
		fill_slots_kernel<<<gridSize, blockSize>>>(slots_d, d_counter, 
			d_alive, rem_parts, rngs_d, opts_d);

		// Retrieve how many particles were actually filled
		int filled = 0;
		cudaMemcpy(&filled, d_counter, sizeof(int), cudaMemcpyDeviceToHost);

		// Retrieve number of alive slots
		cudaMemcpy(&alive_slots, d_alive, sizeof(int), cudaMemcpyDeviceToHost);

		// Update remaining particles
		rem_parts -= filled;
		if (rem_parts < 0) rem_parts = 0;

		// Free memory
		cudaFree(d_counter);
		cudaFree(d_alive);
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
