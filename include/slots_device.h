/**
* @file  slots_device.h
* @brief Contains struct to hold the particle slots for device-side access.
*/
#pragma once

namespace Slots
{
	/**
	* Device-side copy of a Slots object, made by Slots::to_device and freed
	* with free_slots.
	*
	* A plain struct of device pointers, so it can be passed to kernels by
	* value. Each per-particle array has N entries, one per slot, and slot i
	* is handled by thread i in the kernels. Positions and their velocities
	* (x, vx) are in computational coordinates, and vX, vY, vZ are the
	* Cartesian velocity.
	*/
	struct SlotsDevice
	{
		int device_id;		///< ID of the GPU the arrays are on
		int N;				///< Number of slots
		int Z;				///< Impurity atomic number
		double mass;		///< Impurity mass (kg)

		double* pt;			///< Time since the particle was launched (s)
		double* t;			///< Background time of the particle (s)
		double* x;			///< Computational x coordinate
		double* y;			///< Computational y coordinate
		double* z;			///< Computational z coordinate
		int* tidx;			///< Nearest background time index
		int* xidx;			///< x index of the containing grid cell
		int* yidx;			///< y index of the containing grid cell
		int* zidx;			///< z index of the containing grid cell
		double* vx;			///< Computational x velocity (x units/s)
		double* vy;			///< Computational y velocity (y units/s)
		double* vz;			///< Computational z velocity (z units/s)
		double* vX;			///< Cartesian X velocity (m/s)
		double* vY;			///< Cartesian Y velocity (m/s)
		double* vZ;			///< Cartesian Z velocity (m/s)
		int* iters;			///< Number of steps taken. Allocated, but never
							///< initialized or updated on the device
		double* weight;		///< Particle weight, 1 / imp_num at launch
		int* q;				///< Charge state
		int* state;			///< 0 if the particle is alive, 1 if dead

		// Work arrays for the kernels that refill dead slots (see slots.cu)
		bool* all_dead;		///< Single value, true if every slot is dead
		int* dead_indices;	///< Indices of the dead slots, size N
		int* num_dead;		///< Single value, number of entries in
							///< dead_indices
	};
}
