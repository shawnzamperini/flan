/**
* @file  slots_device.h
* @brief Contains struct to hold the particle slots for device-side access.
*/
#pragma once

#include <cstddef>

#include "host_device.h"

namespace Slots
{
	/**
	* Device-side copy of a Slots object, made by Slots::to_device and freed
	* with free_slots.
	*
	* A plain struct of device pointers, so it can be passed to kernels by
	* value. Each per-particle array has N entries, one per slot. The
	* kernels find their slot with slot_of_thread: after fill_slots_gpu
	* has listed the live slots (live_indices), the threads take them in
	* order, per_warp to a warp, and before that, thread k takes slot k.
	* Positions and their velocities
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

		// The plasma at each slot for the collisions, written each step by
		// Collisions::prepare_collisions_gpu (one interpolation for both
		// species) and read by the ion and electron collision kernels.
		// Only set for live, charged slots; coll_nsub_i and coll_nsub_e
		// are 0 for the others.
		double* coll_ne;		///< ne, at least min_ne (m^-3)
		double* coll_te;		///< Te, at least min_te (eV)
		double* coll_ti;		///< Ti, at least min_ti (eV)
		double* coll_uX;		///< Background flow, Cartesian X (m/s)
		double* coll_uY;		///< Background flow, Cartesian Y (m/s)
		double* coll_uZ;		///< Background flow, Cartesian Z (m/s)
		double* coll_bX;		///< B, Cartesian X, where a heat flux is
								///< nonzero (else 0)
		double* coll_bY;		///< B, Cartesian Y (likewise)
		double* coll_bZ;		///< B, Cartesian Z (likewise)
		double* coll_q_par_e;	///< Electron parallel heat flux (W/m^2), 0
								///< without the thermal force
		double* coll_q_par_i;	///< Ion parallel heat flux (W/m^2), likewise
		int* coll_nsub_i;		///< Collision sub-steps with the ions
		int* coll_nsub_e;		///< Collision sub-steps with the electrons

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
		int per_warp;		///< Live slots given to each warp, 1 to 32
							///< (see per_warp_for); 32 while there's
							///< no list
		int resident_warps;	///< Warps of the heavy kernels this GPU runs
							///< at once (set by the first fill)
	};

	/// Threads in a warp
	constexpr int warp_size {32};

	/**
	* How many live slots each warp gets, so that the launch fills the
	* warps the GPU can run at once and no more.
	*
	* The kernels' lanes diverge (most of all in the collisions), and a
	* warp runs its lanes' different paths one after another. While
	* there are fewer live particles than the GPU has resident warps,
	* packing them 32 to a warp would make them take turns in a few warps
	* when they could each have a warp and run at the same time. So they
	* are spread out, as few to a warp as fills the GPU, and packed fully
	* once there are more than it holds, when the warps would wait their
	* turn anyway and packing saves the dead lanes.
	*
	* @param num_live Number of live slots.
	* @param resident_warps Warps the GPU runs at once, 0 if not known.
	* @return Live slots per warp, 1 to warp_size.
	*/
	FLAN_HD inline int per_warp_for(const int num_live,
		const int resident_warps)
	{
		if (resident_warps <= 0) return warp_size;
		const int p {(num_live + resident_warps - 1) / resident_warps};
		return p < 1 ? 1 : (p > warp_size ? warp_size : p);
	}

	/**
	* Position in the live list that thread k takes: the first per_warp
	* lanes of each warp take consecutive entries, and the other lanes
	* none.
	*
	* @param k Thread index over the whole grid.
	* @param per_warp Live slots per warp.
	* @param num_live Number of entries in the list.
	* @return Position in the list, or -1 if the thread has none.
	*/
	FLAN_HD inline int list_position(const int k, const int per_warp,
		const int num_live)
	{
		const int lane {k % warp_size};
		if (lane >= per_warp) return -1;
		const int j {(k / warp_size) * per_warp + lane};
		return j < num_live ? j : -1;
	}

	/**
	* Threads per block for the kernels that use the most registers:
	* ionization and recombination, the collisions and the push (with its
	* sub-steps). At 130 to 255 registers per thread an SM holds only a
	* block or two of 256, so a launch a little over a multiple of the SM
	* count leaves most SMs idle while its last blocks run. Smaller blocks
	* spread the same threads over more SMs; measured, these kernels were
	* as fast or faster with 128. The lighter kernels keep 256.
	*/
	constexpr int heavy_block_size {128};

	/**
	* Number of blocks of block_size threads that cover the slots the
	* kernels run over (num_live of them, per_warp to a warp). Can be 0,
	* which is not a valid launch, so the caller skips the launch then.
	*
	* @param slots_d Particle slots on the device.
	* @param block_size Threads per block, a multiple of warp_size.
	* @return Number of blocks.
	*/
	inline int num_blocks(const SlotsDevice& slots_d, const int block_size)
	{
		const int warps {(slots_d.num_live + slots_d.per_warp - 1)
			/ slots_d.per_warp};
		const int warps_per_block {block_size / warp_size};
		return (warps + warps_per_block - 1) / warps_per_block;
	}

#ifdef __CUDACC__
	/**
	* The slot the calling thread handles, or -1 if it has none.
	*
	* The threads take the live list in order, per_warp to a warp
	* (list_position), so the live particles sit together in warps
	* instead of being scattered among dead slots, as many to a warp as
	* fills the GPU (per_warp_for). Without a list (before the first fill,
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
		const int j {list_position(k, slots_d.per_warp, slots_d.num_live)};
		if (j < 0) return -1;
		return slots_d.live_indices != nullptr ? slots_d.live_indices[j] : j;
	}
#endif
}
