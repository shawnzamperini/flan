#pragma once

#include "background.h"
#include "options.h"
#include "impurity_stats.h"
#include "timer.h"

namespace ImpurityTransport
{
	// To-do: Rename Statistics namespace to just Statistics or something
	Impurity::Statistics follow_impurities(Background::Background& bkg, 
		Options::Options& opts, Timer::Timer& timer);
}
