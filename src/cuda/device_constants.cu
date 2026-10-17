/**
* @file  device_constants.cu
*
* @brief Definitions of the background grid variables kept in GPU constant
*        memory.
*
* The background's time and coordinate arrays are small and read by nearly
* every kernel, so they are kept in constant memory, which is cached and
* faster to read than global memory when the threads in a warp read the same
* address. Background::to_device fills them with cudaMemcpyToSymbol. Constant
* memory is per GPU, and to_device selects the GPU before copying, so each
* GPU gets its own copy.
*
* The arrays have a fixed size of MAX_1D (see device_constants.cuh). Only the
* first tdim, xdim, ydim or zdim elements (dim + 1 for the grid edges) are
* set; the rest are unused.
*
* These are compiled with relocatable device code, so there is one copy of
* each variable per executable. Code that needs them must be linked with the
* same copy that Background::to_device fills in (see tests/CMakeLists.txt).
*/
#include "device_constants.cuh"


// Background times and cell centers. These are the coordinates that fields
// are interpolated between.
__constant__ double d_t[MAX_1D];		///< Background times (s), length tdim
__constant__ double d_x[MAX_1D];		///< x cell centers, length xdim
__constant__ double d_y[MAX_1D];		///< y cell centers, length ydim
__constant__ double d_z[MAX_1D];		///< z cell centers, length zdim

// Cell edges, used to find the cell containing a particle
__constant__ double d_grid_x[MAX_1D];	///< x cell edges, length xdim + 1
__constant__ double d_grid_y[MAX_1D];	///< y cell edges, length ydim + 1
__constant__ double d_grid_z[MAX_1D];	///< z cell edges, length zdim + 1

// Simulation bounds, used by the boundary conditions and when sampling
// starting positions
__constant__ double d_t_min;	///< First background time (s)
__constant__ double d_t_max;	///< Last background time (s)
__constant__ double d_x_min;	///< Lower x bound of the simulation
__constant__ double d_x_max;	///< Upper x bound of the simulation
__constant__ double d_y_min;	///< Lower y bound of the simulation
__constant__ double d_y_max;	///< Upper y bound of the simulation
__constant__ double d_z_min;	///< Lower z bound of the simulation
__constant__ double d_z_max;	///< Upper z bound of the simulation
