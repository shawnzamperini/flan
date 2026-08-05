#pragma once

#include "background.h"
#include "options.h"
#include "impurity_stats.h"
#include "timer.h"

namespace ImpurityTransport
{

	//	
	template <typename T>
	int get_nearest_index_cpu(const std::vector<T>& vec, const T value);

	//
	template <typename T>
	int get_nearest_cell_index_cpu(const std::vector<T>& grid_edges, 
		const T value);

	// To-do: Rename Statistics namespace to just Statistics or something
	Impurity::Statistics follow_impurities(Background::Background& bkg, 
		Options::Options& opts, Timer::Timer& timer);
}
