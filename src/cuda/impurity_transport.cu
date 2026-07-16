#include <cuda_runtime.h>
#include <cstdio>

#include "background_device.h"
#include "device_constants.cuh"
#include "slots_device.h"

namespace ImpurityTransport
{

	// Finding nearest index GPU. This is the same logic as 
	// get_nearest_index_cpu in impurity_transport.cpp just without STL.
	template <typename T>
	__device__ int get_nearest_index_gpu(const T* vec, int N, T value)
	{
		// Manual lower_bound
		int left = 0;
		int right = N;

		// This is essentially std::lower. Find the index of the first value
		// in vec that is larger than value with a binary search.
		while (left < right)
		{
			int mid = (left + right) >> 1;  // (left + right) / 2
			if (vec[mid] < value)
				left = mid + 1;
			else
				right = mid;
		}

		int index = left;

		// Boundary cases
		if (index == 0)
			return 0;

		if (index == N)
			return N - 1;

		// Compare vec[index] vs vec[index-1]
		T diff_low  = vec[index]   - value;
		T diff_prev = vec[index-1] - value;

		// Absolute value
		if (diff_low < 0)  diff_low  = -diff_low;
		if (diff_prev < 0) diff_prev = -diff_prev;

		// Return the closer index
		return (diff_low < diff_prev) ? index : (index - 1);
	}	

	// Finding nearest cell index from a grid on GPU. This is the same logic as 
	// get_nearest_cell_index_cpu in impurity_transport.cpp just without STL.
	template <typename T>
	__device__ int get_nearest_cell_index_gpu(const T* grid_edges, int N, 
		T value)
	{
		// This is essentially std::lower. Find the index of the first value
		// in grid_edges that is larger than value with a binary search.
		int left = 0;
		int right = N;
		while (left < right)
		{
			int mid = (left + right) >> 1;  // (left + right) / 2

			if (grid_edges[mid] < value)
				left = mid + 1;
			else
				right = mid;
		}

		int index = left;

		// Realize that one minus the index represented by lower is the value
		// we're after in the vectors with values at the cell centers.
		//  ____________
		//  |_0_|_1_|_2_|  <-- cell center indices
		//  0   1   2   3  <-- grid_edges indices
		//          ^
		//        lower
		//
		// In this example, we want 1 returned, so we return 2 - 1 = 1. If
		// value is outside the range, return the index of the respective end
		// of the vector.

		// If value is less than everything return the first cell
		if (index == 0)
			return 0;

		// If value is larger than everything return the last cell
		else if (index == N)
			return N - 2;
		else
			return index - 1;
	}

	__global__ void find_containing_cell_kernel(Slots::SlotsDevice slots_d, 
		const Background::BackgroundDevice bkg_d)
	{

		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= slots_d.N) return;

		// d_t in constant memory, included from cuda/device_constants.cuh and
		// lives in Background namespace.
		int tidx {get_nearest_index_gpu(Background::d_t, bkg_d.tdim, 
			slots_d.t[i])};

		// Likewise, d_grid_x is in constant memory
		int xidx {get_nearest_cell_index_gpu(Background::d_grid_x, bkg_d.xdim, 
			slots_d.x[i])};
		int yidx {get_nearest_cell_index_gpu(Background::d_grid_y, bkg_d.ydim, 
			slots_d.y[i])};
		int zidx {get_nearest_cell_index_gpu(Background::d_grid_z, bkg_d.zdim, 
			slots_d.z[i])};

		// Update indices
		slots_d.tidx[i] = tidx;
		slots_d.xidx[i] = xidx;
		slots_d.yidx[i] = yidx;
		slots_d.zidx[i] = zidx;
	}

	// Rewrite in GPU friendly format
	void find_containing_cell_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d)
	{

		// Block and grid size
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;

		// Call kernel to update indices (tidx, xidx, ...) in slots_d
		find_containing_cell_kernel<<<gridSize, blockSize>>>(slots_d, bkg_d);
	}

} // namespace ImpurityTransport
