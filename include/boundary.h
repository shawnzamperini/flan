/**
* @file boundary.h
* @brief Header file for boundary.cpp.
*/
#pragma once

#include <array>
#include <vector>

#include "background.h"
#include "boundary_hd.h"
#include "options.h"
#include "pcg32.h"
#include "slots.h"
#include "zones.h"


// Declared in sputtering.h, which isn't included here to keep this header
// light
namespace Sputtering {class Deposition;}

namespace Boundary
{
	// Fates is in boundary_hd.h, so the GPU code can add to it

	Fates print_fates(const Fates& fates);

	// absorbing_bc, periodic_bc, core_bc, twist_shift_dy, lcfs_bc,
	// core_face and the rest of the per-particle routines are in
	// boundary_hd.h, shared with the GPU

	BoundaryParams make_params(const Background::Background& bkg,
		const Options::Options& opts, const double* chi);

	void check_bounds_cpu(Slots::Slots& slots, 
		const Background::Background& bkg, const Options::Options& opts,
		std::vector<pcg32>& rngs, Fates& fates,
		Sputtering::Deposition* deposition = nullptr);

	void check_bounds_cpu(Slots::Slots& slots, 
		const Background::Background& bkg, const Options::Options& opts,
		std::vector<pcg32>& rngs);

	[[noreturn]] void avalanche_error(int max_gen,
		double chain_multiplication, double mean_yield, double mean_energy);
}
