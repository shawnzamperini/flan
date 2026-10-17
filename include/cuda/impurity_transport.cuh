#pragma once

/**
* @file  impurity_transport.cuh
*
* @brief Header file for impurity_transport.cu
*/
#include <cuda_runtime.h>

#include "background_device.h"
#include "options_device.h"
#include "slots_device.h"


namespace ImpurityTransport
{
	void find_containing_cell_gpu(Slots::SlotsDevice& slots_d,
		const Background::BackgroundDevice& bkg_d);

	void step_gpu(Slots::SlotsDevice& slots_d, const double dt, 
		const Options::OptionsDevice* opts_d, int* max_time_ctr_d);
}
