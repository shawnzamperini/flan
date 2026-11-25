#pragma once
/**
* @file utilities_hd.h
*
* @brief Math helpers shared by the CPU and GPU code (see host_device.h).
*/

#include <cmath>

#include "host_device.h"


namespace Utilities
{
	/**
	* Trilinearly interpolates within a box of eight values.
	*
	* The box has a corner at (x0, y0, z0) and sizes (dx, dy, dz), which may
	* be negative (the GPU stencils point the box toward the neighboring
	* cell, on either side). The value names give the corner's offset,
	* e.g., v100 is at (x0 + dx, y0, z0). The weights are clamped to
	* [0, 1], so a point outside the box gets the value on the nearest face
	* instead of an extrapolation, which can be far off at the edges of a
	* grid. A weight that is NaN (a box of size 0 in some direction) is
	* clamped to 0, so a grid one cell wide in a direction gives that
	* cell's value.
	*
	* Utilities::trilinear_interpolate (the CPU version, taking the far
	* corner) and the GPU stencils both call this, so the two interpolate
	* the same way.
	*
	* @param x0 x of the origin corner.
	* @param y0 y of the origin corner.
	* @param z0 z of the origin corner.
	* @param dx Size of the box in x.
	* @param dy Size of the box in y.
	* @param dz Size of the box in z.
	* @param v000 Value at (x0, y0, z0).
	* @param v100 Value at (x0 + dx, y0, z0).
	* @param v010 Value at (x0, y0 + dy, z0).
	* @param v110 Value at (x0 + dx, y0 + dy, z0).
	* @param v001 Value at (x0, y0, z0 + dz).
	* @param v101 Value at (x0 + dx, y0, z0 + dz).
	* @param v011 Value at (x0, y0 + dy, z0 + dz).
	* @param v111 Value at (x0 + dx, y0 + dy, z0 + dz).
	* @param x x of the point.
	* @param y y of the point.
	* @param z z of the point.
	* @return Interpolated value at (x, y, z).
	*/
	FLAN_HD inline double trilinear_box(
		const double x0, const double y0, const double z0,
		const double dx, const double dy, const double dz,
		const double v000, const double v100, const double v010,
		const double v110, const double v001, const double v101,
		const double v011, const double v111,
		const double x, const double y, const double z)
	{
		// Normalized coordinates, in [0, 1] inside the box. fmax(0, NaN)
		// is 0, so a NaN weight becomes 0.
		const double tx {std::fmin(1.0, std::fmax(0.0, (x - x0) / dx))};
		const double ty {std::fmin(1.0, std::fmax(0.0, (y - y0) / dy))};
		const double tz {std::fmin(1.0, std::fmax(0.0, (z - z0) / dz))};

		// Interpolate along x for the four lower/upper face corners
		const double c00 {v000 + tx * (v100 - v000)};
		const double c01 {v001 + tx * (v101 - v001)};
		const double c10 {v010 + tx * (v110 - v010)};
		const double c11 {v011 + tx * (v111 - v011)};

		// Interpolate along y for the lower and upper edges
		const double c0 {c00 + ty * (c10 - c00)};
		const double c1 {c01 + ty * (c11 - c01)};

		// Interpolate along z
		return c0 + tz * (c1 - c0);
	}

}  // namespace Utilities
