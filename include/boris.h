#pragma once

#include "background.h"
#include "options.h"
#include "slots.h"
#include "slots_device.h"

namespace Boris
{
	void update_velocity_cpu(Slots::Slots& slots, 
		Slots::SlotsDevice& slots_d, const Background::Background& bkg, 
		const Background::BackgroundDevice& bkg_d, const Options::Options& opts,
		const double dt);
}
