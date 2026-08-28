#include <cuda_runtime.h>

#include "background_device.h"
#include "slots_device.h"

namespace Boris
{
	void update_velocity_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d, 
		const double dt);
}
