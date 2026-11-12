/**
* @file impurity_transport.h
* @brief Header file for impurity_transport.cpp.
*/
#pragma once

#include "background.h"
#include "options.h"
#include "impurity_stats.h"
#include "timer.h"

// Declared in sputtering.h, which isn't included here to keep this header
// light
namespace Sputtering {class Source;}

namespace ImpurityTransport
{
	// To-do: Rename Statistics namespace to just Statistics or something
	Impurity::Statistics follow_impurities(Background::Background& bkg, 
		Options::Options& opts, Timer::Timer& timer,
		const Sputtering::Source* source = nullptr);
}
