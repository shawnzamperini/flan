/*
* Launcher for the kernel of test_plasma_interp_cuda: interpolates the
* plasma at points on the GPU with Interpolate::interp_plasma and
* interp_plasma_vec3, the routines the collision, ionization and birth
* kernels use. Kept apart from the test so the test can be compiled by
* the host compiler, which can read the test backgrounds.
*/

#pragma once

#include <vector>

#include "background_device.h"

namespace PlasmaInterpTest
{
	/// A point to interpolate at, with its cell (as in the slots)
	struct Point
	{
		int zone;
		int tidx;
		int xidx;
		int yidx;
		int zidx;
		double t;
		double x;
		double y;
		double z;
	};

	/// Number of values per point: ne, Te, Ti, |E|, uX, uY, uZ
	inline constexpr int num_values {7};

	/**
	* Interpolates the plasma at each point on the GPU.
	*
	* @param bkg_d Background on the device.
	* @param points The points.
	* @return num_values values per point, point after point.
	*/
	std::vector<double> interp_on_gpu(
		const Background::BackgroundDevice& bkg_d,
		const std::vector<Point>& points);
}
