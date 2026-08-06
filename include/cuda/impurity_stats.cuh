#pragma once

#include "slots_device.h"
#include "impurity_stats_device.h"

namespace ImpurityStats
{
	/*
	* @brief CUDA implementation of get_nearest_index
	*/
	template <typename T>
	__device__ int get_nearest_index_gpu(const T* vec, int N, T value);

	/**
	* @brief Wrapper to call record_stats GPU kernel
	*/
	void record_stats_gpu(StatisticsDevice& stats_d, 
		const Slots::SlotsDevice& slots_d, double imp_time_step);
}

