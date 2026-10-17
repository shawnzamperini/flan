#pragma once

/**
* @file  openadas.cuh
*
* @brief Header file for openadas.cu
*/
#include "background_device.h"
#include "openadas_device.h"
#include "pcg32.h"
#include "slots_device.h"


namespace OpenADAS
{
	IonizRecombWarningsDevice alloc_warnings(int device_id);

	void add_warnings_to_host(const IonizRecombWarningsDevice& warnings_d,
		int& ioniz_warnings, int& recomb_warnings, int& max_time_ctr);

	void free_warnings(IonizRecombWarningsDevice& warnings_d);

	void ioniz_recomb_gpu(Slots::SlotsDevice& slots_d,
		const Background::BackgroundDevice& bkg_d,
		const OpenADASDevice& oa_ioniz_d,
		const OpenADASDevice& oa_recomb_d, const double dt,
		const IonizRecombWarningsDevice& warnings_d, pcg32* rngs_d);

}  // OpenADAS
