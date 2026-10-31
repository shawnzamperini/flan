/**
* @file boundary.h
* @brief Header file for boundary.cpp.
*/
#pragma once

#include <array>
#include <vector>

#include "background.h"
#include "options.h"
#include "pcg32.h"
#include "slots.h"
#include "zones.h"


namespace Boundary
{
	/**
	* What happened to particles at the boundaries over a run, counted by
	* check_bounds_cpu. Reaching imp_max_time is counted separately.
	*/
	struct Fates
	{
		long long crossings {};	///< Crossings into another zone
		long long absorbed {};	///< Absorbed at an absorbing face
		long long wall {};		///< Absorbed inside the wall (chi)
		long long core {};		///< Put back in at a core face
		long long time_end {};	///< Reached the end of the background time
		long long runaway {};	///< Absorbed after too many crossings

		Fates& operator+=(const Fates& other);
	};

	void print_fates(const Fates& fates);

	int absorbing_bc(double a, double bound, bool max_bound,
		double buffer = 0.0);

	double periodic_bc(double a, double amin, double amax);

	void core_bc(double& a, double& b, double& c, double a_min,
		double a_max, double a_buffer, double b_min, double b_max,
		double b_buffer, double c_min, double c_max, double c_buffer, 
		bool max_bound, pcg32& rng);

	bool core_face(std::array<double, 3>& pos, int& zone,
		const Zones::Face face, const std::vector<int>& core_zones,
		const Zones::Zones& zones, const std::array<double, 3>& buffers,
		pcg32& rng);

	void check_bounds_cpu(Slots::Slots& slots, 
		const Background::Background& bkg, const Options::Options& opts,
		std::vector<pcg32>& rngs, Fates& fates);

	void check_bounds_cpu(Slots::Slots& slots, 
		const Background::Background& bkg, const Options::Options& opts,
		std::vector<pcg32>& rngs);
}
