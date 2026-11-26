/**
* @file  impurity_stats.cu
*
* @brief CUDA routines for recording the impurity statistics on a GPU.
*
* GPU version of record_stats_cpu in impurity_stats.cpp. Each time step,
* every alive particle adds its contribution to the statistics arrays of
* the cell and time frame it's in. Many particles can be in the same cell,
* so the additions are atomic.
*/
#include <cuda_runtime.h> 
#include <cstdio>

#include "impurity_stats.cuh"
#include "impurity_stats_device.h"
#include "slots_device.h"


namespace ImpurityStats
{

	/**
	* Kernel that adds each alive particle to the statistics of its cell
	* and time frame, one thread per slot.
	*
	* The count is incremented and weight * dt (the occupancy, in number of
	* particles) is added to the weights. The Cartesian and computational
	* velocities and the charge are added multiplied by weight * dt, so they
	* can later be turned into weighted averages. Dead slots are skipped.
	* The cell and time frame are the particle's zone, tidx, xidx, yidx and
	* zidx, so find_containing_cell_gpu must have been called first.
	*
	* @param[in,out] stats_d Impurity statistics on the device. Added to.
	* @param slots_d Particle slots on the device.
	* @param dt Time step (s).
	*/
	__global__ void record_stats_kernel(StatisticsDevice stats_d, 
		Slots::SlotsDevice slots_d, double dt)
	{
		// This thread's slot (the next live one, see slot_of_thread)
		const int i {Slots::slot_of_thread(slots_d)};
		if (i < 0) return;

		// Don't count dead particles
		if (slots_d.state[i]) return;

		// The particle's cell and time frame in its zone, as
		// record_stats_cpu
		const std::size_t idx {index(stats_d, slots_d.zone[i],
			slots_d.tidx[i], slots_d.xidx[i], slots_d.yidx[i],
			slots_d.zidx[i])};

		// Add to stats. atomicAdd is needed since other particles may be
		// adding to the same cell at the same time.
		double p_w_dt {slots_d.weight[i] * dt};  // Change to float
		atomicAdd(&stats_d.counts[idx], 1);
		atomicAdd(&stats_d.weights[idx], p_w_dt);
		atomicAdd(&stats_d.vX[idx], slots_d.vX[i] * p_w_dt);
		atomicAdd(&stats_d.vY[idx], slots_d.vY[i] * p_w_dt);
		atomicAdd(&stats_d.vZ[idx], slots_d.vZ[i] * p_w_dt);
		atomicAdd(&stats_d.vx[idx], slots_d.vx[i] * p_w_dt);
		atomicAdd(&stats_d.vy[idx], slots_d.vy[i] * p_w_dt);
		atomicAdd(&stats_d.vz[idx], slots_d.vz[i] * p_w_dt);
		atomicAdd(&stats_d.q[idx],  slots_d.q[i]  * p_w_dt);
	}
	
	/**
	* Launches record_stats_kernel to add every alive particle on the GPU to
	* the impurity statistics.
	*
	* GPU version of record_stats_cpu. Uses blocks of 256 threads, one
	* thread per slot. The launch is asynchronous, except in DEBUG builds,
	* where it waits for the kernel to finish and prints any error.
	*
	* @param[in,out] stats_d Impurity statistics on the device. Added to.
	* @param slots_d Particle slots on the device, on the same GPU as
	*        stats_d.
	* @param dt Time step (s).
	*/
	void record_stats_gpu(StatisticsDevice& stats_d, 
		const Slots::SlotsDevice& slots_d, double dt)
	{

		// Each stats is assigned to a specific GPU where its data lies
		cudaSetDevice(stats_d.device_id);

		// One thread per live slot (see Slots::slot_of_thread), each
		// adding its particle into stats_d
		int blockSize = 256;
		int gridSize = Slots::num_blocks(slots_d, blockSize);
		if (gridSize == 0) return;	// no live slots
		record_stats_kernel<<<gridSize, blockSize>>>(stats_d, slots_d, dt);

#ifdef DEBUG
		// Check for errors
		cudaError_t err {cudaDeviceSynchronize()};
		if (err != cudaSuccess)
			printf("record_stats_kernel error: %s\n", cudaGetErrorString(err));
#endif

	}  // record_stats_gpu

}  // namespace ImpurityStats
