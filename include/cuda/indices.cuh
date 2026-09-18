namespace Indices
{
	// Finding nearest index CUDA. This is the same logic as 
	// get_nearest_index_cpu in indices.h just without STL.
	template <typename T>
	__device__ int get_nearest_index_cuda(const T* vec, int N, T value)
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


	// Finding nearest cell index from a grid on CUDA. This is the same logic as 
	// get_nearest_cell_index_cpu in indices.h just without STL.
	template <typename T>
	__device__ int get_nearest_cell_index_cuda(const T* grid_edges, int N, 
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

}  // namespace Indices
