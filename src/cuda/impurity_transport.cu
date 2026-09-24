#include <cuda_runtime.h>
#include <cstdio>

#include <curand_kernel.h>

#include "background_device.h"
#include "openadas_device.h"
#include "pcg32.h"
#include "slots_device.h"

#include "device_constants.cuh"
#include "indices.cuh"
#include "utilities.cuh"

namespace ImpurityTransport
{

	__global__ void find_containing_cell_kernel(Slots::SlotsDevice slots_d, 
		const Background::BackgroundDevice bkg_d)
	{

		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;

		// Don't try and access beyond the number of slots (segfault)
		if (i >= slots_d.N) return;

		// If particle isn't alive then skip
		if (slots_d.state[i] > 0) return;

		// d_t in constant memory, included from cuda/device_constants.cuh and
		// lives in Background namespace.
		int tidx {Indices::get_nearest_index_cuda(d_t, bkg_d.tdim, 
			slots_d.t[i])};

		// Likewise, d_grid_x is in constant memory and they have shape dim+1
		int xidx {Indices::get_nearest_cell_index_cuda(d_grid_x, bkg_d.xdim+1, 
			slots_d.x[i])};
		int yidx {Indices::get_nearest_cell_index_cuda(d_grid_y, bkg_d.ydim+1, 
			slots_d.y[i])};
		int zidx {Indices::get_nearest_cell_index_cuda(d_grid_z, bkg_d.zdim+1, 
			slots_d.z[i])};

		// Update indices
		slots_d.tidx[i] = tidx;
		slots_d.xidx[i] = xidx;
		slots_d.yidx[i] = yidx;
		slots_d.zidx[i] = zidx;
	}


	// Find what cell all the particles in slots_d are in, updating the
	// indices accordingly (GPU)
	void find_containing_cell_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d)
	{

		// Block and grid size
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;

		// Call kernel to update indices (tidx, xidx, ...) in slots_d
		find_containing_cell_kernel<<<gridSize, blockSize>>>(slots_d, bkg_d);
	}


	// Step particles
	__global__ void step_kernel(Slots::SlotsDevice slots_d, const double dt)
	{
		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;

		// Don't try and access beyond the number of slots (segfault)
		if (i >= slots_d.N) return;

		// If particle isn't alive then skip
		if (slots_d.state[i] > 0) return;

		slots_d.t[i] += dt;
		slots_d.x[i] += slots_d.vx[i] * dt;
		slots_d.y[i] += slots_d.vy[i] * dt;
		slots_d.z[i] += slots_d.vz[i] * dt;
	}


	// Wrapper to call step_kernel
	void step_gpu(Slots::SlotsDevice& slots_d, const double dt)
	{
		// Block and grid size
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;

		step_kernel<<<gridSize, blockSize>>>(slots_d, dt);

#ifdef DEBUG
		// Check for errors
		cudaError_t err {cudaDeviceSynchronize()};
		if (err != cudaSuccess)
			printf("step_kernel error: %s\n", cudaGetErrorString(err));
#endif
	}

} // namespace ImpurityTransport
