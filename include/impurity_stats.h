#include <vector>

#include "vectors.h"

namespace Impurity
{
	// Class to hold all the vectors related to tracking the statistics of
	// an impurity transport simulation.
	class Statistics
	{
	private:
		Vectors::Vector4D<int> m_counts {};
		Vectors::Vector4D<double> m_weights {};

	public:
		
		// Constructor: Empty
		Statistics(const int dim1, const int dim2, const int dim3, 
			const int dim4)
			: m_counts (Vectors::Vector4D<int> {dim1, dim2, dim3, dim4})
			, m_weights (Vectors::Vector4D<double> {dim1, dim2, dim3, dim4})
		{}

		// Accessors
		Vectors::Vector4D<int> get_counts() {return m_counts;};
		Vectors::Vector4D<double> get_weights() {return m_weights;}
		
	};

}
