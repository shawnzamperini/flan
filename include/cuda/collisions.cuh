#pragma once

/**
* @file  collisions.cuh
*
* @brief Header file for collisions.cu
*/
#include <array>

#include "background_device.h"
#include "collisions_hd.h"
#include "options_device.h"
#include "pcg32.h"
#include "slots_device.h"

#include "impurity_stats.cuh"


namespace Collisions
{
	void collision_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d, const bool elec, 
		const double dt, ImpurityStats::StatisticsDevice& imp_stats_d,
		pcg32* rngs_d, const double elec_mass_amu, const double ion_mass_amu,
		const Options::OptionsDevice* opts_d, unsigned long long* skew_d);

	void prepare_collisions_gpu(Slots::SlotsDevice& slots_d,
		const Background::BackgroundDevice& bkg_d, const double dt,
		const double elec_mass_amu, const double ion_mass_amu,
		const Options::OptionsDevice* opts_d);

	void collisions_gpu(Slots::SlotsDevice& slots_d,
		const Background::BackgroundDevice& bkg_d, const double dt,
		ImpurityStats::StatisticsDevice& imp_stats_d, pcg32* rngs_d,
		const double elec_mass_amu, const double ion_mass_amu,
		const Options::OptionsDevice* opts_d, unsigned long long* skew_d);

	std::array<SkewCounts, 2> skew_counts_from_device(
		const unsigned long long* skew_d, const int device_id);

}  // namespace Collisions
