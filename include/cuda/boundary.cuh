#pragma once

/**
* @file  boundary.cuh
*
* @brief Header file for boundary.cu
*/
#include "background_device.h"
#include "boundary_hd.h"
#include "pcg32.h"
#include "slots_device.h"
#include "sputtering_hd.h"

namespace Boundary
{
	void check_bounds_gpu(Slots::SlotsDevice& slots_d,
		const Background::BackgroundDevice& bkg_d, const BoundaryParams& p,
		pcg32* rngs_d, unsigned long long* fates_d,
		const Sputtering::DepositionDevice* dep_d = nullptr,
		const Sputtering::SelfDevice* self_d = nullptr);

	void add_fates_from_device(const unsigned long long* fates_d,
		const int device_id, Fates& fates);
}
