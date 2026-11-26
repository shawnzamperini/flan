/**
* @file impurity_transport.h
* @brief Header file for impurity_transport.cpp.
*/
#pragma once

#include <cstddef>
#include <string>

#include "background.h"
#include "boundary.h"
#include "options.h"
#include "impurity_stats.h"
#include "timer.h"

// Declared in sputtering.h, which isn't included here to keep this header
// light
namespace Sputtering {class Source; class Deposition;}
namespace Sheath {class Sheath;}

namespace ImpurityTransport
{
	/// GPU memory a run needs on each GPU, by what it's for (bytes)
	struct GpuMemoryNeed
	{
		std::size_t background {};	///< Background::device_bytes
		std::size_t statistics {};	///< Statistics::device_bytes
		std::size_t slots {};		///< Particle slots and their generators
	};

	std::string gpu_memory_shortfall(const GpuMemoryNeed& need,
		std::size_t free_bytes, std::size_t total_bytes, int device_id);

	// To-do: Rename Statistics namespace to just Statistics or something
	Impurity::Statistics follow_impurities(Background::Background& bkg, 
		Options::Options& opts, Timer::Timer& timer,
		const Sputtering::Source* source = nullptr,
		Sputtering::Deposition* deposition = nullptr,
		Boundary::Fates* fates_out = nullptr,
		const Sheath::Sheath* sheath = nullptr);
}
