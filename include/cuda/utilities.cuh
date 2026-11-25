#pragma once

/**
* @file  utilities.cuh
*
* @brief Small math and indexing helpers used in CUDA kernels.
*
* These are the device versions of helpers in utilities.h and vectors.h:
* vector products, neighbor and flattened array indices, and linear and
* bilinear interpolation (trilinear is in utilities_hd.h, shared with the
* CPU). Everything is inline (or a
* template) so it can be used in any kernel without separate compilation.
*/
#include <cmath>
#include <cuda_runtime.h>


namespace Utilities
{

	/**
	* Computes the cross product a x b.
	*
	* out must not be the same array as a or b, since out[0] is written
	* before a and b are done being read.
	*
	* @param a First vector.
	* @param b Second vector.
	* @param[out] out Cross product a x b.
	*/
	__device__ __forceinline__ 
	void cross_product_cuda(const double a[3], const double b[3],
		double out[3])
	{
		out[0] = a[1] * b[2] - a[2] * b[1];
		out[1] = a[2] * b[0] - a[0] * b[2];
		out[2] = a[0] * b[1] - a[1] * b[0];
	}


	/**
	* Computes the dot product of two 3-vectors.
	*
	* @param a First vector.
	* @param b Second vector.
	* @return a . b
	*/
	__device__ __forceinline__ 
	double dot_product_cuda(const double a[3], const double b[3])
	{
		return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	}


	/**
	* Finds the index of the neighboring cell center on val's side of the
	* center of cell idx.
	*
	* This picks the second point to interpolate between. If val is above
	* cell_centers[idx], the neighbor is idx + 1, otherwise (including val
	* exactly at the center) it's idx - 1. At the first or last cell, the
	* only neighbor is returned regardless of side, and the clamped
	* interpolation weights then hold the value constant past the outer
	* cell centers. Same logic as the CPU version, Utilities::
	* get_neighbor_index, written without branches.
	*
	* N must be at least 2. With N = 1 there is no neighbor and the
	* returned index is out of range.
	*
	* @param val Coordinate of the point.
	* @param cell_centers Cell center coordinates, length N.
	* @param idx Index of the cell containing val.
	* @param N Number of cells.
	* @return Index of the neighboring cell.
	*/
	__device__ inline 
	int get_neighbor_index_cuda(const double val,
		const double* cell_centers, int idx, int N)
	{
		// dx > 0 -> +1, dx < 0 -> -1
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

	/**
	* Converts a 4D (t, x, y, z) index to the index into a flattened 4D
	* array, matching Vector4D::calc_index.
	*
	* @param xdim Number of cells in x.
	* @param ydim Number of cells in y.
	* @param zdim Number of cells in z.
	* @param tidx Time index.
	* @param xidx x index.
	* @param yidx y index.
	* @param zidx z index.
	* @return Flattened index.
	*/
	__device__ inline 
	int calc_4d_index_cuda(int xdim, int ydim, 
		int zdim, int tidx, int xidx, int yidx, int zidx)
	{
		// This is from Vector4D::calc_index
		return tidx * (xdim * ydim * zdim) + xidx * (ydim * zdim) 
			+ yidx * zdim + zidx;
	}

	/**
	* Converts a 3D (x, y, z) index to the index into a flattened 3D
	* array, matching Vector3D::calc_index.
	*
	* @param ydim Number of cells in y.
	* @param zdim Number of cells in z.
	* @param xidx x index.
	* @param yidx y index.
	* @param zidx z index.
	* @return Flattened index.
	*/
	__device__ inline 
	int calc_3d_index_cuda(int ydim, int zdim, int xidx, int yidx, int zidx)
	{
		// This is from Vector3D::calc_index
		return xidx * (ydim * zdim) + yidx * zdim + zidx;
	}


	/**
	* Linearly interpolates a fixed-size table at x.
	*
	* Searches the table from the start for the interval containing x, so
	* xarr can be increasing or decreasing. x must be within the range of
	* xarr. If it isn't, yarr[N - 1] is returned, which is only correct if
	* x is past the last entry.
	*
	* @tparam T Floating-point type of the table.
	* @tparam N Length of the table, deduced from the arrays.
	* @param xarr Table x values.
	* @param yarr Table y values.
	* @param x Point to interpolate at.
	* @return Interpolated y value.
	*/
	template <typename T, size_t N>
	__device__
	T linear_interpolate_cuda(const T (&xarr)[N], const T (&yarr)[N],
		T x)
	{
		for (int i = 0; i < N - 1; i++)
		{
			if ((x >= xarr[i] && x <= xarr[i + 1]) ||
				(x <= xarr[i] && x >= xarr[i + 1]))
			{
				const T t = (x - xarr[i]) / (xarr[i + 1] - xarr[i]);
				return yarr[i] + t * (yarr[i + 1] - yarr[i]);
			}
		}

		// Should never happen for valid usage.
		return yarr[N - 1];
	}

	/**
	* Linearly interpolates a fixed-size table with evenly spaced x values
	* at x.
	*
	* Faster than linear_interpolate_cuda for large tables, since the
	* interval containing x is calculated directly instead of searched
	* for. xarr must be increasing and evenly spaced, like the Nanbu table
	* from python/calc_nanbu_A.py. Points outside the table are
	* extrapolated from the first or last interval.
	*
	* @tparam T Floating-point type of the table.
	* @tparam N Length of the table, deduced from the arrays.
	* @param xarr Table x values, increasing and evenly spaced.
	* @param yarr Table y values.
	* @param x Point to interpolate at.
	* @return Interpolated y value.
	*/
	template <typename T, size_t N>
	__device__ inline
	T linear_interpolate_uniform_cuda(const T (&xarr)[N], const T (&yarr)[N],
			T x)
	{
		// Table spacing, from the ends so rounding in the individual
		// entries averages out
		const T dx = (xarr[N - 1] - xarr[0]) / (N - 1);

		// Interval containing x. Clamped so x at (or rounding past) the
		// ends of the table stays inside it.
		int i = static_cast<int>((x - xarr[0]) / dx);
		i = max(0, min(i, static_cast<int>(N) - 2));

		// Use the actual interval width, since the table entries are
		// rounded and not exactly dx apart
		const T t = (x - xarr[i]) / (xarr[i + 1] - xarr[i]);
		return yarr[i] + t * (yarr[i + 1] - yarr[i]);
	}

	/**
	* Performs bilinear interpolation on a flattened two-dimensional
	* lookup table stored in device memory.
	*
	* The function locates the grid cell that contains the point
	* (x, y) with a linear search, computes the interpolation weights in
	* each dimension, and returns the bilinearly interpolated value.
	* Interpolation weights are clamped to the range [0, 1] to prevent
	* extrapolation beyond the boundaries of the tabulated data. nx and ny
	* must both be at least 2.
	*
	* The table is assumed to be stored in row-major order with
	* dimensions nx × ny, where the flattened index corresponding to
	* (ix, iy) is ix * ny + iy.
	*
	* @tparam T Floating-point data type used for the table and coordinate 
	*		  arrays.
	* @param table Flattened two-dimensional lookup table of size nx * ny.
	* @param x_grid Monotonically increasing x-coordinate array of length nx.
	* @param nx Number of x-grid points.
	* @param y_grid Monotonically increasing y-coordinate array of length ny.
	* @param ny Number of y-grid points.
	* @param x X-coordinate at which to evaluate the table.
	* @param y Y-coordinate at which to evaluate the table.
	* @return Bilinearly interpolated value at (x, y).
	*/
	template <typename T>
	__device__ inline T bilinear_interp(const T* __restrict__ table,
		const T* __restrict__ x_grid, int nx, const T* __restrict__ y_grid, 
		int ny, T x, T y)                        
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

		// Prevent extrapolation past edges
		tx = fmin((T)1, fmax((T)0, tx));
		ty = fmin((T)1, fmax((T)0, ty));

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

	// trilinear interpolation is Utilities::trilinear_box in
	// utilities_hd.h, shared with the CPU

}
