/**
* @file boris.h
* @brief Header file for boris.cpp.
*/
#pragma once

#include <array>
#include <vector>

#include "background.h"
#include "options.h"
#include "slots.h"
#include "slots_device.h"

namespace Sheath {class Sheath;}

namespace Boris
{
	void kick(std::array<double, 3>& v, const std::array<double, 3>& E,
		const std::array<double, 3>& B, double q_m, double dt);

	void update_velocity_cpu(Slots::Slots& slots, 
		const Background::Background& bkg, 
		const Options::Options& opts,
		const double dt,
		const Sheath::Sheath* sheath = nullptr,
		const std::vector<int>* skip = nullptr);

	void push_one_cpu(Slots::Slots& slots, const Background::Background& bkg,
		int i, double dt, const Sheath::Sheath* sheath);

	void substep_cpu(Slots::Slots& slots, const Background::Background& bkg,
		double dt, const Sheath::Sheath& sheath, const std::vector<int>& nsub);
}
