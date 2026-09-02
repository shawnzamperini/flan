#pragma once

#include <cuda_runtime.h>

namespace Utilities
{

	// Compute cross-product in-place, storing result in out
	__device__ __forceinline__ 
	void cross_product_cuda(const double a[3], const double b[3],
		double out[3])
	{
		out[0] = a[1] * b[2] - a[2] * b[1];
		out[1] = a[2] * b[0] - a[0] * b[2];
		out[2] = a[0] * b[1] - a[1] * b[0];
	}


	// Compute and return dot product
	__device__ __forceinline__ 
	double dot_product_cuda(const double a[3], const double b[3])
	{
		return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	}


	// Find nearest neighbor index of val using cell_center as the locations
	__device__ inline 
	int get_neighbor_index_cuda(const double val,
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

	// 4D -> 1D index for indexing Vector4Ds that are stored flattened
	__device__ inline 
	int calc_4d_index_cuda(int xdim, int ydim, 
		int zdim, int tidx, int xidx, int yidx, int zidx)
	{
		// This is from Vector4D::calc_index
		return tidx * (xdim * ydim * zdim) + xidx * (ydim * zdim) 
			+ yidx * zdim + zidx;
	}

	// 3D -> 1D index for indexing Vector4Ds that are stored flattened
	__device__ inline 
	int calc_3d_index_cuda(int ydim, int zdim, int xidx, int yidx, int zidx)
	{
		// This is from Vector3D::calc_index
		return xidx * (ydim * zdim) + yidx * zdim + zidx;
	}


	// Bilinearly interpolate the values in flattened 2D arr at (x, y)
	template <typename T>
	__device__ inline T bilinear_interp(
		const T* __restrict__ table,   // flattened 2D array: nx * ny
		const T* __restrict__ x_grid,  // size nx
		int nx,
		const T* __restrict__ y_grid,  // size ny
		int ny,
		T x,                           // physical x value
		T y)                           // physical y value
	{
		// --- 1. Find bracketing x indices ---
		int ix = 0;
		while (ix + 1 < nx && x_grid[ix + 1] < x)
			ix++;

		if (ix >= nx - 1)
			ix = nx - 2;

		T tx = (x - x_grid[ix]) /
					(x_grid[ix + 1] - x_grid[ix]);

		// --- 2. Find bracketing y indices ---
		int iy = 0;
		while (iy + 1 < ny && y_grid[iy + 1] < y)
			iy++;

		if (iy >= ny - 1)
			iy = ny - 2;

		T ty = (y - y_grid[iy]) /
					(y_grid[iy + 1] - y_grid[iy]);

		// --- 3. Compute flattened indices ---
		int idx00 = ix * ny + iy;
		int idx10 = (ix + 1) * ny + iy;
		int idx01 = ix * ny + (iy + 1);
		int idx11 = (ix + 1) * ny + (iy + 1);

		T v00 = table[idx00];
		T v10 = table[idx10];
		T v01 = table[idx01];
		T v11 = table[idx11];

		// --- 4. Bilinear interpolation ---
		T vx0 = v00 + tx * (v10 - v00);
		T vx1 = v01 + tx * (v11 - v01);

		return vx0 + ty * (vx1 - vx0);
	}


	// Trilinearly interpolate the values defined by the 8 vertices at (x, y, z)
	__device__ inline 
	double trilinear_interpolate(
		const double x0, const double y0, const double z0, 
		const double dx, const double dy, const double dz,
		const double v000, const double v100, const double v010, 
		const double v110, const double v001, const double v101, 
		const double v011, const double v111,
		const double x, const double y, const double z)
	{
		// Normalized coordinates in [0,1]
		const double tx = (x - x0) / dx;  // dx = x1 - x0
		const double ty = (y - y0) / dy;
		const double tz = (z - z0) / dz;

		// Interpolate along x for the four lower/upper face corners
		const double c00 = v000 + tx * (v100 - v000);
		const double c01 = v001 + tx * (v101 - v001);
		const double c10 = v010 + tx * (v110 - v010);
		const double c11 = v011 + tx * (v111 - v011);

		// Interpolate along y for the lower and upper edges
		const double c0 = c00 + ty * (c10 - c00);
		const double c1 = c01 + ty * (c11 - c01);

		// Interpolate along z and return
		return c0 + tz * (c1 - c0);
	}

}
