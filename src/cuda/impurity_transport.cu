/**
* @file  impurity_transport.cu
*
* @brief CUDA routines for the main loop steps of impurity_transport.cpp
*        that don't have a file of their own.
*
* GPU versions of find_containing_cell_cpu and step_cpu. The wrappers in
* impurity_transport.cpp call these when use_gpu is on.
*/
#include <cuda_runtime.h>
#include <cstdio>

#include "background_device.h"
#include "openadas_device.h"
#include "options_device.h"
#include "pcg32.h"
#include "slots_device.h"

#include "device_constants.cuh"
#include "indices.cuh"
#include "utilities.cuh"

namespace ImpurityTransport
{

	/**
	* Kernel that finds the background time index and grid cell of each
	* alive particle, one thread per slot.
	*
	* The time index is the nearest background time frame, and the x, y
	* and z indices are the cells containing the particle (see
	* indices.cuh). Dead slots are skipped.
	*
	* @param[in,out] slots_d Particle slots on the device. tidx, xidx, yidx
	*        and zidx are updated.
	* @param bkg_d Background plasma on the device, used for the grid
	*        dimensions.
	*/
	__global__ void find_containing_cell_kernel(Slots::SlotsDevice slots_d, 
		const Background::BackgroundDevice bkg_d)
	{

		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;

		// Don't try and access beyond the number of slots (segfault)
		if (i >= slots_d.N) return;

		// If particle isn't alive then skip
		if (slots_d.state[i] > 0) return;

		// d_t is in constant memory, see device_constants.cuh
		int tidx {Indices::get_nearest_index_cuda(d_t, bkg_d.tdim, 
			slots_d.t[i])};

		// Likewise for the cell edges, which have dim + 1 values
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


	/**
	* Launches find_containing_cell_kernel to update the time and cell
	* indices of every alive particle on the GPU.
	*
	* GPU version of find_containing_cell_cpu. Uses blocks of 256 threads,
	* one thread per slot. The launch is asynchronous, except in DEBUG
	* builds, where it waits for the kernel to finish and prints any error.
	*
	* @param[in,out] slots_d Particle slots on the device. tidx, xidx, yidx
	*        and zidx are updated.
	* @param bkg_d Background plasma on the device.
	*/
	void find_containing_cell_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d)
	{

		// Block and grid size
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;

		// Call kernel to update indices (tidx, xidx, ...) in slots_d
		find_containing_cell_kernel<<<gridSize, blockSize>>>(slots_d, bkg_d);

#ifdef DEBUG
		// Check for errors
		cudaError_t err {cudaDeviceSynchronize()};
		if (err != cudaSuccess)
			printf("find_containing_cell_kernel error: %s\n", 
				cudaGetErrorString(err));
#endif
	}


	/**
	* Kernel that advances each alive particle by one time step, one thread
	* per slot.
	*
	* The time since launch pt and the background time t are increased by
	* dt, and x, y and z are moved by the computational velocity (vx, vy,
	* vz) times dt. If pt is then past the imp_max_time option, the particle
	* is killed and counted in max_time_ctr_d. Dead slots are skipped.
	* Unlike step_cpu, the iteration count isn't updated.
	*
	* @param[in,out] slots_d Particle slots on the device. pt, t, x, y and z
	*        are updated, and state is set to dead for particles past
	*        imp_max_time.
	* @param dt Time step (s).
	* @param opts_d Simulation options on the device, used for
	*        imp_max_time.
	* @param[in,out] max_time_ctr_d Device counter of particles killed for
	*        exceeding imp_max_time. Incremented atomically.
	*/
	__global__ void step_kernel(Slots::SlotsDevice slots_d, const double dt, 
		const Options::OptionsDevice* opts_d, int* max_time_ctr_d)
	{
		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;

		// Don't try and access beyond the number of slots (segfault)
		if (i >= slots_d.N) return;

		// If particle isn't alive then skip
		if (slots_d.state[i] > 0) return;

		slots_d.pt[i] += dt;
		slots_d.t[i] += dt;
		slots_d.x[i] += slots_d.vx[i] * dt;
		slots_d.y[i] += slots_d.vy[i] * dt;
		slots_d.z[i] += slots_d.vz[i] * dt;

		// Check if maximum allowed time has been exceeded, setting to
		// dead if so and incrementing counter
		if (slots_d.pt[i] > opts_d->imp_max_time)
		{
			slots_d.state[i] = 1;
			atomicAdd(max_time_ctr_d, 1);
		}
	}


	/**
	* Launches step_kernel to advance every alive particle on the GPU by one
	* time step.
	*
	* GPU version of step_cpu. Uses blocks of 256 threads, one thread per
	* slot. The launch is asynchronous, except in DEBUG builds, where it
	* waits for the kernel to finish and prints any error.
	*
	* @param[in,out] slots_d Particle slots on the device.
	* @param dt Time step (s).
	* @param opts_d Simulation options on the device, used for
	*        imp_max_time.
	* @param[in,out] max_time_ctr_d Device counter of particles killed for
	*        exceeding imp_max_time, e.g. the max_time member of
	*        OpenADAS::IonizRecombWarningsDevice.
	*/
	void step_gpu(Slots::SlotsDevice& slots_d, const double dt,
		const Options::OptionsDevice* opts_d, int* max_time_ctr_d)
	{
		// Block and grid size
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;

		step_kernel<<<gridSize, blockSize>>>(slots_d, dt, opts_d, 
			max_time_ctr_d);

#ifdef DEBUG
		// Check for errors
		cudaError_t err {cudaDeviceSynchronize()};
		if (err != cudaSuccess)
			printf("step_kernel error: %s\n", cudaGetErrorString(err));
#endif
	}

} // namespace ImpurityTransport
