#include <cuda_runtime.h>

#include "pcg32.h"

#include "slots_device.h"


// Iniitalize RNGs on the device using the tid to select unique seeds
__global__ 
void init_rngs_kernel(pcg32* rngs, uint64_t base_seed) 
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    rngs[tid] = pcg32(base_seed + tid, tid);
}

// Helper function to call init_rngs_kernel
pcg32* init_rngs_cuda(Slots::SlotsDevice& slots_d, uint64_t base_seed)
{

    // This should be moved to an external function for consistency since 
	// every kernel needs the same amount of threads to match the RNG that
	// is assigned to each
	int blockSize = 256;
	int gridSize  = (slots_d.N + blockSize - 1) / blockSize;
	int Nthreads = blockSize * gridSize;

    // Allocate RNG array on device
    pcg32* d_rngs = nullptr;
    cudaMalloc(&d_rngs, Nthreads * sizeof(pcg32));

    // Initialize RNGs on device
    init_rngs_kernel<<<gridSize, blockSize>>>(d_rngs, base_seed);

    // Optional: synchronize for safety during initialization
    cudaDeviceSynchronize();

    return d_rngs;
}
