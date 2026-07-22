#include "device_constants.cuh"

namespace Background
{
	// 1D arrays
    __constant__ double d_t[MAX_1D];
    __constant__ double d_x[MAX_1D];
    __constant__ double d_y[MAX_1D];
    __constant__ double d_z[MAX_1D];
    __constant__ double d_grid_x[MAX_1D];
    __constant__ double d_grid_y[MAX_1D];
    __constant__ double d_grid_z[MAX_1D];

	// 1D array bounds
	__constant__ double d_t_min;
	__constant__ double d_t_max;
	__constant__ double d_x_min;
	__constant__ double d_x_max;
	__constant__ double d_y_min;
	__constant__ double d_y_max;
	__constant__ double d_z_min;
	__constant__ double d_z_max;

	// 3D arrays
    __constant__ double d_dxdX[3 * MAX_1D];
    __constant__ double d_dxdY[3 * MAX_1D];
    __constant__ double d_dxdZ[3 * MAX_1D];
    __constant__ double d_dydX[3 * MAX_1D];
    __constant__ double d_dydY[3 * MAX_1D];
    __constant__ double d_dydZ[3 * MAX_1D];
    __constant__ double d_dzdX[3 * MAX_1D];
    __constant__ double d_dzdY[3 * MAX_1D];
    __constant__ double d_dzdZ[3 * MAX_1D];
}
