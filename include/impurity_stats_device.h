/**
* @file impurity_stats_device.h
* @brief Contains struct to hold impurity statistics data for device-side
*        access.
*/
#pragma once

#include <cstddef>

#include "host_device.h"

namespace ImpurityStats
{
	/**
	* Device-side copy of the impurity statistics arrays, made by
	* Statistics::to_device and freed with free_stats.
	*
	* Each array is the flat data of one of the host's ZonedVectors, which
	* have no ghost cells: zone after zone, each zone (t, x, y, z) on its
	* grid (see index). Each alive particle adds to the cell and time frame
	* it's in every time step (see record_stats_kernel). Apart from counts,
	* every value is added multiplied by the particle's weight * dt, so
	* dividing by weights gives the weighted average.
	*/
	struct StatisticsDevice
	{
		int device_id;		///< ID of the GPU the arrays are on
		int nzones;			///< Number of zones

		/// Start of each zone in the arrays (nzones values, on the device)
		std::size_t* zone_offset;

		int* nx;			///< Cells in x per zone (on the device)
		int* ny;			///< Cells in y per zone (on the device)
		int* nz;			///< Cells in z per zone (on the device)

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

	/**
	* Index into the statistics arrays of a cell and time frame, the same
	* as ZonedVector::calc_index without ghost cells.
	*
	* @param s The statistics on the device.
	* @param zone Zone of the cell.
	* @param t Time frame index.
	* @param i x index of the cell.
	* @param j y index of the cell.
	* @param k z index of the cell.
	* @return Index into the arrays.
	*/
	FLAN_HD inline std::size_t index(const StatisticsDevice& s,
		const int zone, const int t, const int i, const int j, const int k)
	{
		const std::size_t ny {static_cast<std::size_t>(s.ny[zone])};
		const std::size_t nz {static_cast<std::size_t>(s.nz[zone])};
		const std::size_t frame {static_cast<std::size_t>(s.nx[zone])
			* ny * nz};
		return s.zone_offset[zone] + static_cast<std::size_t>(t) * frame
			+ static_cast<std::size_t>(i) * ny * nz
			+ static_cast<std::size_t>(j) * nz + static_cast<std::size_t>(k);
	}
}
