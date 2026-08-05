#pragma once

#include "background_device.h"
#include "slots_device.h"

namespace Boundary
{
	/**
	* @brief Check time and spatial boundaries and handle accordingly (GPU)
	*/
	void check_bounds_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d, 
		const int tbound_type_int, 
		const int min_xbound_type_int, const int max_xbound_type_int, 
		const int min_ybound_type_int, const int max_ybound_type_int, 
		const int min_zbound_type_int, const int max_zbound_type_int, 
		const double imp_xbound_buffer, 
		const double imp_ybound_buffer, 
		const double imp_zbound_buffer, 
		const double lcfs_x);
}
