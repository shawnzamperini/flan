#pragma once

#include "background.h"
#include "options.h"
#include "impurity_stats.h"
#include "timer.h"

namespace ImpurityTransport
{
	// Structure used in intializing particles. We don't actually use individual
	// particle logic outside of this to keep it SoA.
	struct ParticleInit 
	{
		double x, y, z;
		double vx, vy, vz;
		double vX, vY, vZ;
	};

	// To-do: Rename Statistics namespace to just Statistics or something
	Impurity::Statistics follow_impurities(Background::Background& bkg, 
		Options::Options& opts, Timer::Timer& timer);
}
