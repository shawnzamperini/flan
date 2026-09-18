#include <algorithm>
#include <cmath>
#include <vector>


namespace Indices
{
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


	template <typename T>
	int get_nearest_cell_index_cpu(const std::vector<T>& grid_edges, 
		const T value)
	{
		// Get the index of the first value in grid_edges that is larger
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
		// iterator that points past the last element of a vector, so to return
		// the cell we need to subtract by 2. 
		else if (lower == grid_edges.end()) 
		{
			return index - 2;
		}

		//else return lower - grid_edges.begin() - 1;
		else return index - 1;
	}

}  // namespace Indices
