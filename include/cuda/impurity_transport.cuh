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
	* @brief Check time and spatial boundaries and handle accordingly (GPU)
	*/
	void check_bounds_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d, 
		const int tbound_type_int, const double imp_xbound_buffer, 
		const int min_xbound_type_int, const double lcfs_x);

}
