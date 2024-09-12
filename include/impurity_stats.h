#ifndef IMPURITY_STATS_H
#define IMPURITY_STATS_H

#include <vector>

#include "vectors.h"

namespace Impurity
{
	// Class to hold all the vectors related to tracking the statistics of
	// an impurity transport simulation.
	// counts: The number of times an impurity was registered as being in a
	//   given cell at each timestep.
	// weights: The accumulated particle weight in each cell. 
	class Statistics
	{
	private:
		int m_dim1 {};
		int m_dim2 {};
		int m_dim3 {};
		int m_dim4 {};
		int m_size {};
		Vectors::Vector4D<int> m_counts {};
		Vectors::Vector4D<double> m_weights {};

	public:
		
		// Constructor: Empty
		Statistics(const int dim1, const int dim2, const int dim3, 
			const int dim4)
			: m_dim1 {dim1}
			, m_dim2 {dim2}
			, m_dim3 {dim3}
			, m_dim4 {dim4}
			, m_counts (Vectors::Vector4D<int> {dim1, dim2, dim3, dim4})
			, m_weights (Vectors::Vector4D<double> {dim1, dim2, dim3, dim4})
		{}
		
		/*
		Statistics(const Statistics&) 
		{
			std::cout << "Copy constructor\n";
		}

		Statistics(Statistics&&) noexcept 
		{
			std::cout << "Move constructor\n";
		}
		
		Statistics& operator=(const Statistics&) 
		{
			std::cout << "Copy assignment\n";
			return *this;
		}

		Statistics& operator=(Statistics&&) noexcept 
		{
			std::cout << "Move assignment\n";
			return *this;
		}
		*/

		// Accessors
		Vectors::Vector4D<int>& get_counts() {return m_counts;};
		Vectors::Vector4D<double>& get_weights() {return m_weights;}

		// Functions to increase counts
		void add_counts(const int tidx, const int xidx, const int yidx, 
			const int zidx, const int value = 1)
		{
			m_counts(tidx, xidx, yidx, zidx) += value;
		}

		// Fucntion to increase weights
		void add_weights(const int tidx, const int xidx, const int yidx, 
			const int zidx, const int value = 1.0)
		{
			m_weights(tidx, xidx, yidx, zidx) += value;
		}

		// Adding two Statistics together is okay. It is defined here as just
		// the counts and weights together.
		Statistics operator+ (const Statistics& other) const
		{
			Statistics ret_stats {m_dim1, m_dim2, m_dim3, m_dim4};
			ret_stats.m_counts = m_counts + other.m_counts;
			ret_stats.m_weights = m_weights + other.m_weights;

			return ret_stats;
		}

		
	};

}

#endif
