#pragma once

#include <cuda_runtime.h>

namespace Utilities
{

	// Compute cross-product in-place, storing result in out
	__device__ inline void cross_product_cuda(const double a[3], const double b[3],
		double out[3])
	{
		out[0] = a[1] * b[2] - a[2] * b[1];
		out[1] = a[2] * b[0] - a[0] * b[2];
		out[2] = a[0] * b[1] - a[1] * b[0];
	}


	// Compute and return dot product
	__device__ inline double dot_product_cuda(const double a[3], const double b[3])
	{
		return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	}


	// Find nearest neighbor index of val using cell_center as the locations
	__device__ inline int get_neighbor_index_cuda(const double val,
		const double* cell_centers, int idx, int N)
	{
		// dx > 0 → +1, dx < 0 → -1
		double dx_from_center = val - cell_centers[idx];
		int side = 2 * (dx_from_center > 0.0) - 1;

		// Check grid edges
		int is_left_edge  = (idx == 0);
		int is_right_edge = (idx == N - 1);

		// Offset logic identical to CPU version
		int offset =
			  side * (1 - is_left_edge - is_right_edge)
			+ 1    * is_left_edge
			+ (-1) * is_right_edge;

		return idx + offset;
	}
}
