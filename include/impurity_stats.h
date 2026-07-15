/**
* @file impurity_stats.h
*
* @brief Header file for impurity_stats.cpp
*/

#pragma once

#include <numeric>
#include <vector>

#include "background.h"
#include "flan_types.h"
#include "impurity.h"
#include "impurity_stats_device.h"
#include "options.h"
#include "slots.h"
#include "vectors.h"

#ifdef USE_CUDA
#include <cuda_runtime.h> 
#endif

namespace Impurity
{

	void free_stats_d(ImpurityStats::StatisticsDevice& stats_d, int device_id);

	class Statistics
	{
	private:
		int m_dim1 {};
		int m_dim2 {};
		int m_dim3 {};
		int m_dim4 {};
		int m_size {};
		std::vector<float> m_track_t {};
		std::vector<float> m_track_x {};
		std::vector<float> m_track_y {};
		std::vector<float> m_track_z {};
		std::vector<float> m_track_vx {};
		std::vector<float> m_track_vy {};
		std::vector<float> m_track_vz {};
		Vectors::Vector4D<int> m_counts {};
		Vectors::Vector4D<BkgFPType> m_weights {};
		Vectors::Vector4D<BkgFPType> m_density {};
		Vectors::Vector4D<BkgFPType> m_vX {};
		Vectors::Vector4D<BkgFPType> m_vY {};
		Vectors::Vector4D<BkgFPType> m_vZ {};
		Vectors::Vector4D<BkgFPType> m_vx {};
		Vectors::Vector4D<BkgFPType> m_vy {};
		Vectors::Vector4D<BkgFPType> m_vz {};
		Vectors::Vector4D<BkgFPType> m_charge {};
		Vectors::Vector4D<BkgFPType> m_s {};

	public:
		
		// Constructor
		Statistics();
		Statistics(const int dim1, const int dim2, const int dim3, 
			const int dim4);

		// Accessors
		const int get_dim1() const;
		const int get_dim2() const;
		const int get_dim3() const;
		const int get_dim4() const;
		std::vector<float>& get_track_t();
		std::vector<float>& get_track_x();
		std::vector<float>& get_track_y();
		std::vector<float>& get_track_z();
		std::vector<float>& get_track_vx();
		std::vector<float>& get_track_vy();
		std::vector<float>& get_track_vz();
		Vectors::Vector4D<int>& get_counts();
		Vectors::Vector4D<BkgFPType>& get_weights();
		Vectors::Vector4D<BkgFPType>& get_density();
		Vectors::Vector4D<BkgFPType>& get_vX();
		Vectors::Vector4D<BkgFPType>& get_vY();
		Vectors::Vector4D<BkgFPType>& get_vZ();
		Vectors::Vector4D<BkgFPType>& get_vx();
		Vectors::Vector4D<BkgFPType>& get_vy();
		Vectors::Vector4D<BkgFPType>& get_vz();
		Vectors::Vector4D<BkgFPType>& get_charge();
		Vectors::Vector4D<BkgFPType>& get_s();

		// Overload of + to add counts and weights together, returned as a 
		// new Statistics object.
		Statistics operator+ (const Statistics& other) const;

		// Functions to increase counts
		void add_counts(const int tidx, const int xidx, const int yidx, 
			const int zidx, const int value);

		// Function to increase weights
		void add_weights(const int tidx, const int xidx, const int yidx, 
			const int zidx, const BkgFPType value);

		// Function add each velocity component to the corresponding array
		// location
		void add_vels(const int tidx, const int xidx, const int yidx,
			const int zidx, const BkgFPType vX, const BkgFPType vY, 
			const BkgFPType vZ, const BkgFPType vx, const BkgFPType vy,
			const BkgFPType vz, const Background::Background& bkg);

		// Function to increase charge
		void add_charge(const int tidx, const int xidx, const int yidx, 
			const int zidx, const BkgFPType value);

		// Function to increase s
		void add_s(const int tidx, const int xidx, const int yidx, 
			const int zidx, const BkgFPType value);

		// Update particle track.
		void update_track(const Impurity& imp);
 
		// Pack all the arrays into a single buffer that can be sent across
		// MPI processes
		void pack(std::vector<double>& buf) const;

		// Unpack all the arrays that were packed via pack
		void unpack(const std::vector<double>& buf);

		// Calculate the impurity density using the data stored in counts and 
		// weights.
		void calc_density(const Background::Background& bkg, 
			const int tot_imp_num, const double imp_source_scale_fact = 1.0);

		// Calculate the average velocity in each cell. 
		void calc_vels();

		// Calculate the average charge. 
		void calc_charge();

		// Calculate the average Nanbu collisionality strength. 
		void calc_s();

		// GPU transfer
		ImpurityStats::StatisticsDevice to_device(int device_id);

		// Reduce GPU statistics into host-side clas
		void add_stats_device(ImpurityStats::StatisticsDevice& stats_d, 
			int device_id);
	};

	// Helper function to reduce a Statistic object across MPI ranks
	Statistics reduce_stats(const Statistics& local_stats);

	// Free GPU memory
	void free_stats(ImpurityStats::StatisticsDevice& stats_d, int device_id);

	// Record stats using particle values in slots.
	void record_stats_cpu(Statistics& imp_stats, const Slots::Slots& slots, 
		const Options::Options& opts, const double imp_time_step);

}
