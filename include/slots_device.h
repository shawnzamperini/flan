/**
* @file  slots_device.h
* @brief Contains struct to hold the particle slots for device-side access.
*/
#pragma once

#include <cstddef>

namespace Slots
{
	/**
	* Device-side copy of a Slots object, made by Slots::to_device and freed
	* with free_slots.
	*
	* A plain struct of device pointers, so it can be passed to kernels by
	* value. Each per-particle array has N entries, one per slot. The
	* kernels find their slot with slot_of_thread: after fill_slots_gpu
	* has listed the live slots (live_indices), thread k handles the k-th
	* live slot, and before that, slot k. Positions and their velocities
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
		int* zone;			///< Zone the particle is in. The coordinates
							///< and cell indices are in this zone.
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
		double* prompt_by;	///< Particle time up to which a return to the
							///< wall is prompt, below 0 for none (see
							///< Slots::prompt_by)
		int* generation;	///< Times self-sputtered (Slots::generation)
		int* nsub;			///< Near-wall sub-steps of each slot this step
							///< (set when the sheath is on, see
							///< Boris::update_velocity_gpu)

		// Work arrays for the kernels that refill dead slots (see slots.cu)
		bool* all_dead;		///< Single value, true if every slot is dead
		int* dead_indices;	///< Indices of the dead slots, size N
		int* num_dead;		///< Single value, number of entries in
							///< dead_indices
		void* select_temp;	///< Work space for cub::DeviceSelect, which
							///< lists the dead slots (allocated by
							///< fill_slots_gpu when first needed)
		std::size_t select_temp_bytes;	///< Size of select_temp in bytes

		// The live slots, listed by fill_slots_gpu after filling, so the
		// kernels run one thread per live slot rather than per slot
		int* live_indices;	///< Indices of the live slots in increasing
							///< order, size N, or nullptr before the
							///< first fill (then every slot is used)
		int* num_live_d;	///< Single value, number of entries in
							///< live_indices (written by the selection)
		int num_live;		///< The same on the host, which sizes the
							///< launches; N while there's no list
	};

	/**
	* Number of blocks of block_size threads that cover the slots the
	* kernels run over (num_live of them). Can be 0, which is not a valid
	* launch, so the caller skips the launch then.
	*
	* @param slots_d Particle slots on the device.
	* @param block_size Threads per block.
	* @return Number of blocks.
	*/
	inline int num_blocks(const SlotsDevice& slots_d, const int block_size)
	{
		return (slots_d.num_live + block_size - 1) / block_size;
	}

#ifdef __CUDACC__
	/**
	* The slot the calling thread handles, or -1 if it has none.
	*
	* Thread k (over the whole grid) takes the k-th entry of the live
	* list, so the live particles fill whole warps instead of being
	* scattered among dead slots. Without a list (before the first fill,
	* or in tests that set up slots by hand) thread k takes slot k. A slot
	* in the list can still be dead: particles that die during a step stay
	* listed until the next fill, so the kernels keep checking the state.
	*
	* @param slots_d Particle slots on the device.
	* @return Slot index, or -1 past the last slot.
	*/
	__device__ inline int slot_of_thread(const SlotsDevice& slots_d)
	{
		const int k {static_cast<int>(blockIdx.x * blockDim.x
			+ threadIdx.x)};
		if (k >= slots_d.num_live) return -1;
		return slots_d.live_indices != nullptr ? slots_d.live_indices[k] : k;
	}
#endif
}
