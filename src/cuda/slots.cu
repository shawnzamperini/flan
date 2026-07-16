#include "slots.cuh"
#include "slots_device.h"

#include <cstdio>


namespace Slots
{
	
	// Initialize a new particle and return it
	__device__ ParticleInitDevice make_new_particle_device()
	{
		ParticleInitDevice p;

		// To-do
		p.x = 0.0;
		p.y = 0.0;
		p.z = 0.0;
		p.vx = 0.0;
		p.vy = 0.0;
		p.vz = 0.0;
		p.vX = 0.0;
		p.vY = 0.0;
		p.vZ = 0.0;
		p.weight = 1.0;
		p.q = 0.0;

		return p;
	}

	/**
	* @brief GPU kernel to fill slots, replacing dead particles with alive ones.
	*/
	//__global__ void fill_slots_kernel(double* x, double* y, double* z,
	//	double* vx, double* vy, double* vz, double* vX, double* vY, double* vZ,
	//	double* weight, int* q, int* state, int N, int* counter,int rem_parts)
	__global__ void fill_slots_kernel(SlotsDevice slots_d, int* counter, 
		int rem_parts)
	{
		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= slots_d.N) return;

		// Skip alive slots
		if (slots_d.state[i] == 0) return;

		// Atomically claim a particle index. This prevents the race condition
		// in which multiple threads grab the same index.
		int idx = atomicAdd(counter, 1);

		// If we ran out of particles, stop
		if (idx >= rem_parts) return;

		// Create a new particle (device-side initializer)
		ParticleInitDevice p = make_new_particle_device();

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
	}

	// Replace dead particles with alive ones, as long as remaining particles
	// are greater than zero.
	void fill_slots_gpu(SlotsDevice& slots_d, int& rem_parts)
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

		// Call GPU kernel to fill slots.
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;

		//fill_slots_kernel<<<gridSize, blockSize>>>(
		//	slots_d.x, slots_d.y, slots_d.z,
		//	slots_d.vx, slots_d.vy, slots_d.vz, slots_d.vX, slots_d.vY, 
		//	slots_d.vZ, slots_d.weight, slots_d.q, slots_d.state, slots_d.N, 
		//	d_counter, rem_parts);
		fill_slots_kernel<<<gridSize, blockSize>>>(slots_d, d_counter, 
			rem_parts);

		// Retrieve how many particles were actually filled
		int filled = 0;
		cudaMemcpy(&filled, d_counter, sizeof(int), cudaMemcpyDeviceToHost);

		// Update remaining particles
		rem_parts -= filled;
		if (rem_parts < 0) rem_parts = 0;

		// Free memory
		cudaFree(d_counter);
	}
}
