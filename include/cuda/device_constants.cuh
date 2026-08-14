#pragma once


// Values from each allocated block below.
// Total memory:  = 27.3 KB

// As of the time of this writing, CUDA GPUs have 64 KB of constant memory
// (this is apparently true for all of them). These numbers are set with 
// that in mind. If you change, please update the comments below that 
// calculate the memory allocated.
constexpr int MAX_1D {200};

// Memory: 7 * MAX_1D * double = 7 * 100 * 8 B = 5,600 B = 5.6 KB
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
