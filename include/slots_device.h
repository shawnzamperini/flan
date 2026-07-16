#pragma once

namespace Slots
{
    // POD struct for GPU-side data access
    struct SlotsDevice
    {
		int device_id;
        int N;
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
    };
}
