#include <cuda_runtime.h>

#include "background_device.h"
#include "slots_device.h"

namespace ImpurityTransport
{

	/**
	*
	*/
	template <typename T>
	__device__ int get_nearest_index_gpu(const T* vec, int N, T value);

	/**
	*
	*/
	template <typename T>
	__device__ int get_nearest_cell_index_hd(const T* grid_edges, int N, 
		T value);
	
	/**
	*
	*/
	void find_containing_cell_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d);

	/**
	* @brief Wrapper to call init_rng_kernel, which initializes the RNG for 
	*   each slot on the device
	*/
	void init_slot_rngs(Slots::SlotsDevice& slots_d, 
		const unsigned long long seed);

}
