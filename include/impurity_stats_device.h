/**
* @file impurity_stats_device.h
* @brief Contains struct to hold impurity statistics data for device-side
*        access.
*/
#pragma once

namespace ImpurityStats
{
	/**
	* Device-side copy of the impurity statistics arrays, made by
	* Statistics::to_device and freed with free_stats.
	*
	* Each array is 4D, (t, x, y, z) on the background grid, flattened the
	* same way as the background (see calc_4d_index_cuda). Each alive
	* particle adds to the cell and time frame it's in every time step (see
	* record_stats_kernel). Apart from counts, every value is added
	* multiplied by the particle's weight * dt, so dividing by weights gives
	* the weighted average.
	*/
	struct StatisticsDevice
	{
		int device_id;		///< ID of the GPU the arrays are on
		int Nt;				///< Number of background times
		int Nx;				///< Number of cells in x
		int Ny;				///< Number of cells in y
		int Nz;				///< Number of cells in z

		int *counts;		///< Number of particle time steps recorded
		double *weights;	///< Sum of weight * dt (occupancy)

		double *vX;			///< Sum of Cartesian X velocity * weight * dt
		double *vY;			///< Sum of Cartesian Y velocity * weight * dt
		double *vZ;			///< Sum of Cartesian Z velocity * weight * dt

		double *vx;			///< Sum of computational x velocity * weight
							///< * dt
		double *vy;			///< Sum of computational y velocity * weight
							///< * dt
		double *vz;			///< Sum of computational z velocity * weight
							///< * dt

		double *s;			///< Sum of the Nanbu collision parameter s *
							///< weight * dt (added in collisions.cu)

		double *q;			///< Sum of charge state * weight * dt
	};
}
