#pragma once

/**
* @file  indices.cuh
*
* @brief CUDA routines for finding the index of the time frame or grid cell
*        that a value falls in.
*
* These are the device versions of the routines in indices.h, with the same
* results. std::lower_bound can't be used in device code, so the binary
* search is written out by hand. They are templates, so they are defined
* here in the header.
*/

namespace Indices
{
	/**
	* Returns the index of the value in a sorted array that is nearest to
	* value.
	*
	* Used to find the nearest background time frame. Same as
	* get_nearest_index_cpu in indices.h. If value is exactly halfway
	* between two values, the lower index is returned. Values outside the
	* range of vec return the index of the first or last element.
	*
	* @tparam T Floating-point data type of the array and value.
	* @param vec Sorted array of values, e.g., the background times (d_t).
	* @param N Length of vec.
	* @param value Value to find the nearest index for.
	* @return Index of the nearest value in vec.
	*/
	template <typename T>
	__device__ int get_nearest_index_cuda(const T* vec, int N, T value)
	{
		// Manual lower_bound
		int left = 0;
		int right = N;

		// Same as std::lower_bound. Binary search for the index of the first
		// value in vec that is not less than value.
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

		// Compare vec[index] (the first value not less than value) and
		// vec[index-1] (the last one less than value)
		T diff_low  = vec[index]   - value;
		T diff_prev = vec[index-1] - value;

		// Absolute value
		if (diff_low < 0)  diff_low  = -diff_low;
		if (diff_prev < 0) diff_prev = -diff_prev;

		// Return the closer index, or the lower one if they're equally close
		return (diff_low < diff_prev) ? index : (index - 1);
	}	


	/**
	* Returns the index of the grid cell containing value.
	*
	* Same as get_nearest_cell_index_cpu in indices.h. A value exactly on
	* the edge between two cells is placed in the lower cell. Values
	* outside the grid return the index of the first or last cell.
	*
	* @tparam T Floating-point data type of the grid edges and value.
	* @param grid_edges Sorted grid edges, e.g., d_grid_x.
	* @param N Number of grid edges, one more than the number of cells
	*        (e.g., bkg_d.xdim + 1).
	* @param value Value to find the containing cell for.
	* @return Index of the cell containing value, between 0 and N - 2.
	*/
	template <typename T>
	__device__ int get_nearest_cell_index_cuda(const T* grid_edges, int N, 
		T value)
	{
		// Same as std::lower_bound. Binary search for the index of the first
		// value in grid_edges that is not less than value.
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

		// If value is larger than everything return the last cell. There
		// are N - 1 cells, so the last one is N - 2.
		else if (index == N)
			return N - 2;
		else
			return index - 1;
	}

}  // namespace Indices
