#pragma once

/**
* @file  slots.cuh
*
* @brief Header file for slots.cu
*/
#include "slots_device.h"

namespace Slots
{

	/**
	* A newly created particle, returned by make_new_particle_cuda and then
	* written into a slot.
	*
	* t is in s, x, y and z are computational coordinates, vx, vy and vz
	* are computational velocities, vX, vY and vZ are Cartesian velocities
	* (m/s), weight is the particle weight, the idx members are the time and
	* cell indices and q is the charge state.
	*/
	struct ParticleInitDevice {
		double t, x, y, z;
		double vx, vy, vz;
		double vX, vY, vZ;
		double weight;
		int tidx, xidx, yidx, zidx, q;
	};

	void fill_slots_gpu(SlotsDevice& slots_d, int& rem_parts, int& alive_slots,
		pcg32* rngs_d, Options::OptionsDevice* opts_d, 
		const Background::BackgroundDevice& bkg_d);

	bool all_dead_gpu(SlotsDevice& slots_d);
	
}
