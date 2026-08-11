#pragma once

#include "slots_device.h"

namespace Slots
{

	// Structure to hold initialized particle
	struct ParticleInitDevice {
		double t, x, y, z;
		double vx, vy, vz;
		double vX, vY, vZ;
		double weight;
		int tidx, xidx, yidx, zidx, q;
	};

	/**
	* @brief Replace dead particles with alive ones, as long as remaining 
	*   particles are greater than zero (GPU).
	*
	* @param slots Slots object to swap out dead particles in
	* @param rem_parts Number of remaining alive particles to track
	*/
	void fill_slots_gpu(SlotsDevice& slots_d, int& rem_parts);

	/**
	* @brief Check if slots contain all dead particles
	*
	* @param slots_d Device-side slots struct to check
	*
	* @return true if all dead, false if not
	*/
	bool all_dead_gpu(SlotsDevice& slots_d);
	
}
