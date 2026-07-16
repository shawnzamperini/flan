#include <cuda_runtime.h> 
#include <cstdio>

#include "impurity_stats.cuh"
#include "impurity_stats_device.h"
#include "slots_device.h"


namespace ImpurityStats
{

	__global__ void record_stats_kernel(StatisticsDevice stats_d, 
		Slots::SlotsDevice slots_d, double imp_time_step)
	{
		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;

		// Avoid out of bounds indexing
		if (i >= slots_d.N) return;

		// 4D to 1D index
		//int Nt {stats_d.Nt};
		int Nx {stats_d.Nx};
		int Ny {stats_d.Ny};
		int Nz {stats_d.Nz};
		int idx {slots_d.tidx[i] * (Nx * Ny * Nz) + slots_d.xidx[i] 
			* (Ny * Nz) + slots_d.yidx[i] * Nz + slots_d.zidx[i]};

		// Add to stats. Need to use atomicAdd here to avoid a race condition :(
		double p_w {slots_d.weight[i]};
		atomicAdd(&stats_d.counts[idx], 1);
		atomicAdd(&stats_d.weights[idx], p_w * imp_time_step);
		atomicAdd(&stats_d.vX[idx], slots_d.vX[i] * p_w);
		atomicAdd(&stats_d.vY[idx], slots_d.vY[i] * p_w);
		atomicAdd(&stats_d.vZ[idx], slots_d.vZ[i] * p_w);
		atomicAdd(&stats_d.vx[idx], slots_d.vx[i] * p_w);
		atomicAdd(&stats_d.vy[idx], slots_d.vy[i] * p_w);
		atomicAdd(&stats_d.vz[idx], slots_d.vz[i] * p_w);
		atomicAdd(&stats_d.q[idx],  slots_d.q[i]  * p_w);
	}
	
	void record_stats_gpu(StatisticsDevice& stats_d, 
		const Slots::SlotsDevice& slots_d, double imp_time_step)
	{

		// Each stats is assigned to a specific GPU where its data lies
		cudaSetDevice(stats_d.device_id);

		// Call GPU kernel to record stats. We use the number of slots to
		// calculate the grid size because we effectively want to loop through
		// all the slots and aggregate their data into stats_d "one slot at 
		// a time"
		int blockSize = 256;
		int gridSize = (slots_d.N + blockSize - 1) / blockSize;
		record_stats_kernel<<<gridSize, blockSize>>>(stats_d, slots_d, 
			imp_time_step);

		// Check for errors
		cudaError_t err {cudaDeviceSynchronize()};
		if (err != cudaSuccess)
			printf("record_stats_kernel error: %s\n", cudaGetErrorString(err));

	}

}
