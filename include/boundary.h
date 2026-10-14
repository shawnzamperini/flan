/**
* @file boundary.h
* @brief Header file for boundary.cpp.
*/
#pragma once

#include <vector>

#include "background.h"
#include "options.h"
#include "pcg32.h"
#include "slots.h"


namespace Boundary
{

	// Test for absorbing boundary
	int absorbing_bc(double a, double bound, bool max_bound,
		double buffer = 0.0);

	// Test for periodic boundary
	double periodic_bc(double a, double amin, double amax);

	// Test for core boundary
	void core_bc(double& a, double& b, double& c, double a_min,
		double a_max, double a_buffer, double b_min, double b_max,
		double b_buffer, double c_min, double c_max, double c_buffer, 
		bool max_bound, pcg32& rng);
		
	// Top-level function for checking all boundary conditions
	void check_bounds_cpu(Slots::Slots& slots, 
		const Background::Background& bkg, const Options::Options& opts,
		std::vector<pcg32>& rngs);
}
