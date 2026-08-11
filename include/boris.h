#pragma once

#include "background.h"
#include "options.h"
#include "slots.h"
#include "slots_device.h"

namespace Boris
{
	void update_velocity_cpu(Slots::Slots& slots, 
		const Background::Background& bkg, 
		const Options::Options& opts,
		const double dt);
}
