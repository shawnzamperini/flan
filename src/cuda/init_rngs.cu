/**
* @file  init_rngs.cu
*
* @brief CUDA routines for creating the random number generators used on a
*        GPU.
*
* Each slot gets its own PCG32 random number generator, stored in a device
* array indexed by slot. Kernels that need random numbers (collisions,
* ionization/recombination, boundary conditions and filling slots) use the
* generator of the slot they are working on.
*/
#include <cuda_runtime.h>

#include "pcg32.h"
#include "slots_device.h"


/**
* Kernel that creates one PCG32 random number generator per thread.
*
* Thread tid gets seed base_seed + tid and stream tid, so every generator
* in the array produces a different sequence.
*
* @param rngs Device array to fill. Must hold one generator per thread in
*        the launch.
* @param base_seed Base seed shared by every generator in the array.
*/
__global__
void init_rngs_kernel(pcg32* rngs, uint64_t base_seed) 
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    rngs[tid] = pcg32(base_seed + tid, tid);
}

/**
* Allocates and initializes the random number generators for the slots on
* a GPU.
*
* Launches init_rngs_kernel with blocks of 256 threads, so the array is
* rounded up to a whole number of blocks and has at least slots_d.N
* generators. Generator i belongs to slot i. Waits for the kernel to
* finish before returning. The array is allocated on the current device,
* so that must be the GPU slots_d is on.
*
* @param slots_d Particle slots on the device, used for the number of
*        slots.
* @param base_seed Base seed for the generators. See init_rngs_kernel.
* @return Device pointer to the array of generators. The caller owns it
*         and must free it with free_rngs.
*/
pcg32* init_rngs_cuda(Slots::SlotsDevice& slots_d, uint64_t base_seed)
{
	// Select device these slots are on
	cudaSetDevice(slots_d.device_id);

	int blockSize = 256;
	int gridSize  = (slots_d.N + blockSize - 1) / blockSize;
	int Nthreads = blockSize * gridSize;

    // Allocate RNG array on device
    pcg32* d_rngs = nullptr;
    cudaMalloc(&d_rngs, Nthreads * sizeof(pcg32));

    // Initialize RNGs on device
    init_rngs_kernel<<<gridSize, blockSize>>>(d_rngs, base_seed);

    // Wait for the generators to be ready before returning
    cudaDeviceSynchronize();

    return d_rngs;
}


/**
* Frees the random number generators allocated by init_rngs_cuda.
*
* @param rngs_d Device array of generators. Set to nullptr.
* @param device_id ID of the GPU the generators are on.
*/
void free_rngs(pcg32*& rngs_d, int device_id)
{
	cudaSetDevice(device_id);
	cudaFree(rngs_d);
	rngs_d = nullptr;
}
