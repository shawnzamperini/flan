/**
* @file collisions.h
* @brief Header file for collisions.cpp
*/
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <tuple>
#include <vector>

#include "background.h"
#include "constants.h"
#include "impurity_stats.h"
#include "options.h"
#include "pcg32.h"
#include "slots.h"
#include "collisions_hd.h"

namespace Collisions
{
	void nanbu_coll(Slots::Slots& slots, const Background::Background& bkg,
		const Options::Options& opts, bool elec, const double dt,
		Impurity::Statistics& imp_stats, std::vector<pcg32>& rngs);

	double nanbu_s_coeff(const int q, const double T, const double Te,
		const double ne, const double imp_mass, const double mass_kg);

	/**
	* @name Sub-cycling of the Nanbu collisions
	*
	* Each call, Nanbu's model collides the impurity with one background
	* particle and turns their relative velocity g by an angle chi, with
	* <cos chi> = exp(-s) and s proportional to dt / g^3. The friction and
	* the thermal force on a slow, heavy impurity come mostly from slow
	* background particles, since those have the largest s.
	*
	* That is also where a long time step goes wrong. Once s >> 1 for a
	* pair, the deflection is fully random and a longer step can't make
	* that one collision count for more. In reality the impurity would have
	* met many slow background particles in that time, each pushing on it.
	* So the forces per unit time come out weaker than they should, more
	* so the longer the step.
	*
	* How long the step is for the collisions is measured by s0, the value
	* of s for a background particle at its thermal speed v_th = sqrt(T /
	* m) (see nanbu_s0). For another particle s = s0 / w^3, with w its
	* relative speed in units of v_th, so everything slower than s0^(1/3)
	* thermal speeds has s > 1. Relative to their values for small steps
	* (the fluid results), the forces on a slow, heavy impurity are about:
	*
	*   s0      ion thermal force   electron thermal force   friction
	*   0.01    0.95                0.98                     0.97
	*   0.03    0.91                0.96                     0.94
	*   0.1     0.83                0.91                     0.88
	*   1       0.51                0.67                     0.63
	*   10      0.11                0.20                     0.25
	*
	* The shortfall shrinks slowly, as about s0^(2/3). s0 goes as
	* ne q^2 dt / T^(3/2), and for electrons it is about sqrt(m_i / m_e) ~
	* 60 times what it is for ions, so cold, dense plasma, high charge
	* states and electrons are where it matters.
	*
	* Shrinking the time step would fix it but slow everything else down.
	* Instead, only the collision is split: n sub-steps of dt / n, each with
	* a newly sampled background particle, so the impurity meets n of them.
	* Each sub-step then has s0 / n. The plasma at the impurity location is
	* interpolated once for all of them.
	*/
	///@{

	/**
	* s0: the Nanbu s term for a background particle moving at its thermal
	* speed v_th = sqrt(T / m) relative to the impurity. It measures how
	* collisional a time step is for an impurity that is slow compared
	* with the background particles; a fast impurity has larger relative
	* speeds and so smaller s, which makes s0 the careful choice.
	*
	* @param s_coeff From nanbu_s_coeff.
	* @param dt Time step (s).
	* @param T Temperature of the background species (eV).
	* @param mass_kg Mass of the background species (kg).
	* @return s0 (dimensionless).
	*/
	inline double nanbu_s0(const double s_coeff, const double dt,
		const double T, const double mass_kg)
	{
		const double v_th {std::sqrt(T * Constants::ev_to_j / mass_kg)};
		return s_coeff * dt / (v_th * v_th * v_th);
	}

	/**
	* Number of sub-steps to split one Nanbu collision into, so that each
	* has s0 at or below target: s0 / target rounded up, and at most cap.
	*
	* The target sets how close the forces get to their small-step values
	* (see the table above: 0.03 gives 0.91 or better). The cap limits the
	* cost in very collisional plasma, where s0 can be in the hundreds;
	* past s0 = target * cap the forces fall off again, but an impurity
	* there is already tied closely to the background flow.
	*
	* If s0 is at or below target, or isn't a number, this is 1: the plain
	* model, with the same random numbers as without sub-cycling. cap = 1
	* turns sub-cycling off.
	*
	* @param s0 From nanbu_s0, for the full time step.
	* @param target Largest s0 a sub-step is meant to have
	*        (imp_coll_subcycle_target).
	* @param cap Largest number of sub-steps (imp_coll_subcycle_cap).
	* @return Number of sub-steps, from 1 to cap.
	*/
	inline int nanbu_num_substeps(const double s0, const double target,
		const int cap)
	{
		// Also catches NaN, for which the comparison is false
		if (!(s0 > target)) return 1;

		// Compare as doubles first: s0 / target can be far too large for
		// an int in very collisional plasma
		const double n {std::ceil(s0 / target)};
		if (n >= cap) return cap;
		return std::max(static_cast<int>(n), 1);
	}
	///@}


	std::tuple<double, double, double> sample_bkg_velocity(
		const double T, const double uX, const double uY, const double uZ,
		const double m, const double eps, const SkewShape& shape,
		const double bX, const double bY, const double bZ, pcg32& rng,
		SkewCounts& counts);

	SkewCounts get_skew_counts(const bool elec);

	void reset_skew_counts();

	void add_skew_counts(const std::array<SkewCounts, 2>& counts);

	std::array<SkewCounts, 2> reduce_skew_counts();

	void check_thermal_force(const Background::Background& bkg,
		const Options::Options& opts);

	void print_skew_counts(const std::array<SkewCounts, 2>& totals);

	// The heat-flux skew of the background distribution (thermal force):
	// SkewShape, skew_shape, skew_eps, skew_factor, sample_bkg, etc. are
	// in collisions_hd.h, shared with the GPU.
}
