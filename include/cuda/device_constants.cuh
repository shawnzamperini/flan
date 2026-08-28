#pragma once


// As of the time of this writing, CUDA GPUs have 64 KB of constant memory
// (this is apparently true for all of them). These numbers are set with 
// that in mind. If you change, please update the comments below that 
// calculate the memory allocated.
constexpr int MAX_1D {500};

// Small 1D arrays
extern __constant__ double d_t[MAX_1D];
extern __constant__ double d_x[MAX_1D];
extern __constant__ double d_y[MAX_1D];
extern __constant__ double d_z[MAX_1D];
extern __constant__ double d_grid_x[MAX_1D];
extern __constant__ double d_grid_y[MAX_1D];
extern __constant__ double d_grid_z[MAX_1D];

// 1D array bounds
extern __constant__ double d_t_min;
extern __constant__ double d_t_max;
extern __constant__ double d_x_min;
extern __constant__ double d_x_max;
extern __constant__ double d_y_min;
extern __constant__ double d_y_max;
extern __constant__ double d_z_min;
extern __constant__ double d_z_max;
