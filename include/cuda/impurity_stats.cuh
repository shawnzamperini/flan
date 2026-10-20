#pragma once

/**
* @file  impurity_stats.cuh
*
* @brief Header file for impurity_stats.cu
*/
#include "slots_device.h"
#include "impurity_stats_device.h"

namespace ImpurityStats
{
	void record_stats_gpu(StatisticsDevice& stats_d, 
		const Slots::SlotsDevice& slots_d, double dt);
}
