/**
* @file  indices.h
*
* @brief Routines for finding the index of the time frame or grid cell that a
*        value falls in.
*
* These are templates, so they are defined here in the header instead of in a
* .cpp file. The CUDA versions are in indices.cuh.
*/
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>


namespace Indices
{
	/**
	* Returns the index of the value in a sorted vector that is nearest to
	* value.
	*
	* Used to find the nearest background time frame. If value is exactly
	* halfway between two values, the lower index is returned. Values
	* outside the range of vec return the index of the first or last
	* element.
	*
	* @tparam T Floating-point data type of the vector and value.
	* @param vec Sorted vector of values, e.g., the background times.
	* @param value Value to find the nearest index for.
	* @return Index of the nearest value in vec.
	*/
	template <typename T>
	int get_nearest_index_cpu(const std::vector<T>& vec, const T value)
	{
		auto lower = std::lower_bound(vec.begin(), vec.end(), value);

		if (lower == vec.begin())
		{
			return 0;
		}
		if (lower == vec.end())
		{
			return vec.size() - 1;
		}

		auto prev = lower - 1;
		if (std::fabs(*lower - value) < std::fabs(*prev - value))
		{
			return lower - vec.begin();
		}
		else
		{
			return prev - vec.begin();
		}
	}

	/**
	* Returns the index of the grid cell containing value.
	*
	* A value exactly on the edge between two cells is placed in the lower
	* cell. Values outside the grid return the index of the first or last
	* cell.
	*
	* @tparam T Floating-point data type of the grid edges and value.
	* @param grid_edges Sorted grid edges, one more than the number of cells.
	* @param value Value to find the containing cell for.
	* @return Index of the cell containing value, between 0 and
	*         grid_edges.size() - 2.
	*/
	template <typename T>
	int get_nearest_cell_index_cpu(const std::vector<T>& grid_edges,
		const T value)
	{
		// Get the index of the first value in grid_edges that is not smaller
		// than value.
		auto lower = std::lower_bound(grid_edges.begin(), grid_edges.end(),
			value);

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
		int index = std::distance(grid_edges.begin(), lower);

		// If less than everything, return the first cell.
		if (lower == grid_edges.begin()) return 0;

		// If larger than everything, return the last cell. Note end() is the
		// iterator that points past the last element of a vector, so to
		// return the cell we need to subtract by 2.
		else if (lower == grid_edges.end())
		{
			return index - 2;
		}

		//else return lower - grid_edges.begin() - 1;
		else return index - 1;
	}
}  // namespace Indices
