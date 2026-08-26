#pragma once

#ifdef USE_CUDA
#include <curand_kernel.h>
#endif

namespace Slots
{
    // POD struct for GPU-side data access
    struct SlotsDevice
    {
		int device_id;
        int N;
        int Z;
		double mass;
        double* t;
        double* x;
        double* y;
        double* z;
        int* tidx;
        int* xidx;
        int* yidx;
        int* zidx;
        double* vx;
        double* vy;
        double* vz;
        double* vX;
        double* vY;
        double* vZ;
        double* weight;
        int* q;
        int* state;
		bool* all_dead;

#ifdef USE_CUDA
		
		// Each slot gets it RNG, otherwise we could get non-independent
		// random number streams
		curandState* rng;
#endif

    };
}
