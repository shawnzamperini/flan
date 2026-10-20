#pragma once

/**
* @file  device_constants.cuh
*
* @brief Header file for device_constants.cu
*/


/**
* Size of each array in constant memory.
*
* CUDA GPUs have 64 KB of constant memory. With seven arrays of MAX_1D
* doubles plus eight scalars, MAX_1D = 500 uses about 28 KB. The grid edge
* arrays hold dim + 1 values, so each background dimension must be less than
* MAX_1D.
*/
constexpr int MAX_1D {500};

// Background times and cell centers
extern __constant__ double d_t[MAX_1D];
extern __constant__ double d_x[MAX_1D];
extern __constant__ double d_y[MAX_1D];
extern __constant__ double d_z[MAX_1D];

// Cell edges
extern __constant__ double d_grid_x[MAX_1D];
extern __constant__ double d_grid_y[MAX_1D];
extern __constant__ double d_grid_z[MAX_1D];

// Simulation bounds
extern __constant__ double d_t_min;
extern __constant__ double d_t_max;
extern __constant__ double d_x_min;
extern __constant__ double d_x_max;
extern __constant__ double d_y_min;
extern __constant__ double d_y_max;
extern __constant__ double d_z_min;
extern __constant__ double d_z_max;
