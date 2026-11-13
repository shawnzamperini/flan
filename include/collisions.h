/**
* @file collisions.h
* @brief Header file for collisions.cpp
*/
#pragma once

#include <cmath>
#include <tuple>
#include <vector>

#include "background.h"
#include "constants.h"
#include "impurity_stats.h"
#include "options.h"
#include "pcg32.h"
#include "slots.h"

namespace Collisions
{
	void nanbu_coll(Slots::Slots& slots, const Background::Background& bkg,
		const Options::Options& opts, bool elec, const double dt,
		Impurity::Statistics& imp_stats, std::vector<pcg32>& rngs);

	/**
	* Counts from sampling the skewed (Grad) ion distribution, summed over
	* every collision on this MPI rank since the start of the run (or the
	* last reset_skew_counts). They show how well the skew was reproduced:
	* draws where the Grad factor is negative were dropped (clipped), so
	* negative / draws is the share of the Maxwellian that was clipped.
	*/
	struct SkewCounts
	{
		long long samples {};		///< Background ions sampled with a skew
		long long draws {};			///< Maxwellian draws, accepted or not
		long long negative {};		///< Draws where the Grad factor was < 0
		long long above_bound {};	///< Draws where it was above the
									///< rejection bound (|w| > 5)
		long long eps_capped {};	///< Samples where |eps| > 1 was capped
	};

	double skew_bound(const double eps);

	std::tuple<double, double, double> sample_bkg_velocity(
		const double T, const double uX, const double uY, const double uZ,
		const double m, const double eps, const double bX, const double bY,
		const double bZ, pcg32& rng, SkewCounts& counts);

	SkewCounts get_skew_counts();

	void reset_skew_counts();

	/**
	* @name Heat-flux skew of the background distribution (thermal force)
	*
	* A Maxwellian is symmetric about the flow velocity. With a parallel
	* heat flux q_par, it isn't: in the frame moving with the flow, the
	* fast particles lean toward the direction q_par points (they carry the
	* heat) and the slow ones lean the other way (so there's no net flow).
	* Since Coulomb collisions are strongest with slow particles (the rate
	* goes as 1 / g^3), an impurity colliding with this distribution gets a
	* net push opposite to q_par, i.e., toward higher temperature: the
	* thermal force. Sampling the colliding background particle from the
	* skewed distribution gives it with no separate force term.
	*
	* The skew is the lowest-order (Grad 13-moment) correction to the
	* Maxwellian,
	*
	*   f = f_M [1 + eps w_par (w^2 / 5 - 1)],
	*
	* where w = (v - u) / v_th is the velocity relative to the flow in units
	* of the thermal speed v_th = sqrt(T / m), w_par its component along B,
	* and eps the size of the skew (grad_skew_eps). It has the same density,
	* flow and temperature as f_M, and its heat flux is exactly q_par.
	*/
	///@{

	/**
	* Size of the heat-flux skew, eps = q_par / (n T v_th), v_th =
	* sqrt(T / m).
	*
	* eps is the heat flux relative to n T v_th, roughly the largest heat
	* flux the plasma could carry (its thermal energy moving at the thermal
	* speed). eps = 0 is a Maxwellian, and |eps| ~ 0.1 is a mild skew. Its
	* sign is q_par's: positive when heat flows along B.
	*
	* @param q_par Parallel conductive heat flux (W/m^2), positive along B.
	* @param n Density of the species (m^-3).
	* @param T Temperature of the species (eV).
	* @param mass_kg Mass of the species (kg).
	* @return eps (dimensionless).
	*/
	inline double grad_skew_eps(const double q_par, const double n,
		const double T, const double mass_kg)
	{
		const double T_J {T * Constants::ev_to_j};
		const double v_th {std::sqrt(T_J / mass_kg)};
		return q_par / (n * T_J * v_th);
	}

	/**
	* Factor the Grad distribution differs from the Maxwellian by at a
	* velocity, 1 + eps w_par (w^2 / 5 - 1) (see grad_skew_eps). Above
	* w^2 = 5 (fast particles) it favors moving along the heat flux, below
	* it (slow particles) against it. It is negative for large enough |eps|
	* w^3, i.e., far in the tails or for strong skews, where the expansion
	* breaks down.
	*
	* @param eps Size of the skew, from grad_skew_eps.
	* @param w_par Velocity relative to the flow along B, in units of
	*        v_th.
	* @param w_sq Square of the velocity relative to the flow, in units of
	*        v_th^2.
	* @return The factor f / f_M.
	*/
	inline double grad_factor(const double eps, const double w_par,
		const double w_sq)
	{
		return 1.0 + eps * w_par * (w_sq / 5.0 - 1.0);
	}
	///@}
}
