#pragma once

#include "background_device.h"
#include "pcg32.h"
#include "slots_device.h"

namespace Boundary
{

	// Function to check for absorbing boundary condition
	__device__ 
	int absorbing_bc_gpu(double a, double bound, bool max_bound,
		double buffer = 0.0);

	// Function to check for perioidic boundary condition
	__device__ 
	double periodic_bc_gpu(double a, double amin, double amax);

	// Function to check for core boundary condition
	__device__
	void core_bc_gpu(double& a, double& b, double& c, double a_min,
		double a_max, double a_buffer, double b_min, double b_max, 
		double b_buffer, double c_min, double c_max, double c_buffer, 
		bool max_bound, pcg32& rng);

	/**
	* @brief Check time and spatial boundaries and handle accordingly (GPU)
	*/
	void check_bounds_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d, 
		const int tbound_type_int, 
		const int min_xbound_type_int, const int max_xbound_type_int, 
		const int min_ybound_type_int, const int max_ybound_type_int, 
		const int min_zbound_type_int, const int max_zbound_type_int, 
		const double imp_xbound_buffer, 
		const double imp_ybound_buffer, 
		const double imp_zbound_buffer, 
		const double lcfs_x, pcg32* rngs_d);
}
