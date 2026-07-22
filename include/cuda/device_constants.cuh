#pragma once


// Values from each allocated block below.
// Total memory: 2.73 + 10.55 = 13.28 KB

namespace Background 
{
	// As of the time of this writing, CUDA GPUs have 64 KB of constant memory
	// (this is apparently true for all of them). These numbers are set with 
	// that in mind. If you change, please update the comments below that 
	// calculate the memory allocated.
	constexpr int MAX_1D {50};

	// Memory: 7 * MAX_1D * double = 7 * 50 * 8 B = 1,600 B = 2.73 KB
    extern __constant__ double d_t[MAX_1D];
    extern __constant__ double d_x[MAX_1D];
    extern __constant__ double d_y[MAX_1D];
    extern __constant__ double d_z[MAX_1D];
    extern __constant__ double d_grid_x[MAX_1D];
    extern __constant__ double d_grid_y[MAX_1D];
    extern __constant__ double d_grid_z[MAX_1D];

	// Memory _ * double = __ B 
	// 1D array bounds
	extern __constant__ double d_t_min;
	extern __constant__ double d_t_max;
	extern __constant__ double d_x_min;
	extern __constant__ double d_x_max;
	extern __constant__ double d_y_min;
	extern __constant__ double d_y_max;
	extern __constant__ double d_z_min;
	extern __constant__ double d_z_max;

	// Memory: 9 * (3 * MAX_1D) * double = 9 * 150 * 8 B = 10.5 KB
	// The basis vectors are 3D arrays, so we need to allocate 3 * MAX_1D bytes. 
    extern __constant__ double d_dxdX[3 * MAX_1D];
    extern __constant__ double d_dxdY[3 * MAX_1D];
    extern __constant__ double d_dxdZ[3 * MAX_1D];
    extern __constant__ double d_dydX[3 * MAX_1D];
    extern __constant__ double d_dydY[3 * MAX_1D];
    extern __constant__ double d_dydZ[3 * MAX_1D];
    extern __constant__ double d_dzdX[3 * MAX_1D];
    extern __constant__ double d_dzdY[3 * MAX_1D];
    extern __constant__ double d_dzdZ[3 * MAX_1D];
}
