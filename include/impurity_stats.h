/**
* @file impurity_stats.h
*
* @brief Header file for impurity_stats.cpp
*/

#ifndef IMPURITY_STATS_H
#define IMPURITY_STATS_H

#include <numeric>
#include <vector>

#include "vectors.h"
#include "impurity.h"
#include "background.h"

namespace Impurity
{
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
		Vectors::Vector4D<double> m_density {};
		Vectors::Vector4D<double> m_vx {};
		Vectors::Vector4D<double> m_vy {};
		Vectors::Vector4D<double> m_vz {};
		Vectors::Vector4D<double> m_gyrorad {};
		bool m_vel_stats {};

	public:
		
		// Constructor
		Statistics(const int dim1, const int dim2, const int dim3, 
			const int dim4, const bool vel_stats);
		
		// Accessors
		Vectors::Vector4D<int>& get_counts();
		Vectors::Vector4D<double>& get_weights();
		Vectors::Vector4D<double>& get_density();
		Vectors::Vector4D<double>& get_vx();
		Vectors::Vector4D<double>& get_vy();
		Vectors::Vector4D<double>& get_vz();
		Vectors::Vector4D<double>& get_gyrorad();
		bool get_vel_stats();

		// Overload of + to add counts and weights together, returned as a 
		// new Statistics object.
		Statistics operator+ (const Statistics& other) const;

		// Functions to increase counts
		void add_counts(const int tidx, const int xidx, const int yidx, 
			const int zidx, const int value = 1);

		// Function to increase weights
		void add_weights(const int tidx, const int xidx, const int yidx, 
			const int zidx, const int value = 1.0);

		// Function add each velocity component to the corresponding array
		// location
		void add_vels(const int tidx, const int xidx, const int yidx,
			const int zidx, const double vx, const double vy, const double vz);

		// Function add the calculated gyroradius to the corresponding array
		// location
		void add_gyrorad(const int tidx, const int xidx, 
			const int yidx, const int zidx, const Impurity& imp, 
			const Background::Background& bkg);

		// Calculate the impurity density using the data stored in counts and 
		// weights.
		void calc_density(const Background::Background& bkg, 
			const int tot_imp_num, const double imp_source_scale_fact = 1.0);

		// Calculate the average velocity in each cell. 
		void calc_vels();

		// Calculate the average gyroradius. 
		void calc_gyrorad();
	};

}

#endif
