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

	/**
	* Counts from sampling a skewed background distribution, summed
	* over every collision with that species on this MPI rank since the
	* start of the run (or the last reset_skew_counts). They show how well the skew was reproduced:
	* draws where the skew factor is negative were dropped (clipped), so
	* negative / draws is the share of the Maxwellian that was clipped.
	*/
	struct SkewCounts
	{
		long long samples {};		///< Background particles sampled with a skew
		long long draws {};			///< Maxwellian draws, accepted or not
		long long negative {};		///< Draws where the skew factor was < 0
		long long above_bound {};	///< Draws where it was above the
									///< rejection bound (|w| > 5)
		long long eps_capped {};	///< Samples where |eps| > 1 was capped
	};

	struct SkewShape;

	std::tuple<double, double, double> sample_bkg_velocity(
		const double T, const double uX, const double uY, const double uZ,
		const double m, const double eps, const SkewShape& shape,
		const double bX, const double bY, const double bZ, pcg32& rng,
		SkewCounts& counts);

	SkewCounts get_skew_counts(const bool elec);

	void reset_skew_counts();

	std::array<SkewCounts, 2> reduce_skew_counts();

	void check_thermal_force(const Background::Background& bkg,
		const Options::Options& opts);

	void print_skew_counts(const std::array<SkewCounts, 2>& totals);

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
	* The skewed distribution is
	*
	*   f = f_M [1 + eps w_par (a + b w^2 + c w^4)],
	*
	* where w = (v - u) / v_th is the velocity relative to the flow in units
	* of the thermal speed v_th = sqrt(T / m), w_par its component along B,
	* eps the size of the skew (skew_eps) and a, b, c its shape (SkewShape):
	* at which speeds the particles lean along the heat flux and at which
	* against it. Any shape used here has the same density, flow and
	* temperature as f_M and a heat flux of exactly q_par. What the shape
	* changes is how hard the skew pushes an impurity: the same heat flux
	* carried by a lean among faster particles pushes less, since fast
	* particles collide less.
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
	inline double skew_eps(const double q_par, const double n,
		const double T, const double mass_kg)
	{
		const double T_J {T * Constants::ev_to_j};
		const double v_th {std::sqrt(T_J / mass_kg)};
		return q_par / (n * T_J * v_th);
	}

	/**
	* Shape of the skew: the coefficients of f / f_M = 1 + eps w_par
	* (a + b w^2 + c w^4), and bound, the largest |w_par (a + b w^2 +
	* c w^4)| for speeds up to 5 v_th, used for rejection sampling (see
	* skew_bound). Made by skew_shape.
	*/
	struct SkewShape
	{
		double a {};
		double b {};
		double c {};
		double bound {};
	};

	/**
	* Shape of the skew with a given push on impurities.
	*
	* The three coefficients are fixed by three conditions on the skewed
	* distribution (averages over the Maxwellian, per unit eps):
	*   1. No flow, <w_par^2 (a + b w^2 + c w^4)> = 0: the flow is already
	*      the Maxwellian's drift.
	*   2. The heat flux is q_par, 1/2 <w^2 w_par^2 (...)> = 1.
	*   3. The push on a slow, heavy impurity, which is proportional to the
	*      "friction moment" <w_par^2 (...) / w^3>, is push times what the
	*      lowest-order (Grad 13-moment) shape w_par (w^2 / 5 - 1) gives.
	* They give c = (1 - push) / 25, b = 1/5 - 14 c and a = -1 + 35 c.
	* push = 1 is the Grad shape (c = 0). push < 1 moves the lean to faster
	* particles, which collide less.
	*
	* push is a property of the background's distribution, not of the
	* impurity. For a fluid background it comes from the closure (see
	* skew_push_ion and skew_push_elec). For a kinetic background it could be
	* computed from the distribution itself.
	*
	* @param push Push on a slow, heavy impurity relative to the Grad
	*        shape's, at the same heat flux.
	* @return The shape.
	*/
	inline SkewShape skew_shape(const double push)
	{
		SkewShape s {};
		s.c = (1.0 - push) / 25.0;
		s.b = 0.2 - 14.0 * s.c;
		s.a = -1.0 + 35.0 * s.c;

		// Largest |w (a + b w^2 + c w^4)| for w up to 5 (w_par = w is the
		// worst direction), on a fine grid
		constexpr double w_max {5.0};
		constexpr int n {2000};
		for (int i {}; i <= n; ++i)
		{
			const double w {w_max * i / n};
			const double w2 {w * w};
			s.bound = std::max(s.bound,
				std::abs(w * (s.a + s.b * w2 + s.c * w2 * w2)));
		}
		return s;
	}

	/**
	* Push of the ion skew, for a collisional (fluid) background.
	*
	* With the heat flux q = -kappa grad T and Braginskii's ion conductivity
	* kappa = 3.906 n T tau_i / m_i, the Grad shape pushes a slow, heavy
	* impurity of charge Z with a thermal force of 0.6 sqrt(2) 3.906 Z^2
	* grad T = 3.31 Z^2 grad T. The accurate result in that limit, from
	* the full distortion of the ion distribution, is 2.65 Z^2 grad T
	* (Chapman 1958; the limit of the coefficient beta_i in Stangeby's
	* book). The Grad shape puts too much of the lean on slow ions.
	*/
	inline constexpr double skew_push_ion {2.6517
		/ (0.6 * 1.4142135623730951 * 3.906)};

	/**
	* Push of the electron skew, for a collisional (fluid) background.
	*
	* With Braginskii's electron conductivity kappa = 3.16 n T tau_e / m_e,
	* the Grad shape gives a thermal force of 0.6 * 3.16 Z^2 grad T = 1.90
	* Z^2 grad T, but Braginskii's is 0.71 Z^2 grad T: collisions with ions
	* slow fast electrons much less than slow ones, so the real
	* (Spitzer-Harm) distortion sits at high speeds. The values 3.16 and
	* 0.71 are for a hydrogenic plasma (Z_eff = 1).
	*/
	inline constexpr double skew_push_elec {0.71 / (0.6 * 3.16)};

	/**
	* Factor the skewed distribution differs from the Maxwellian by at a
	* velocity, 1 + eps w_par (a + b w^2 + c w^4). Fast particles favor
	* moving along the heat flux, and slow ones against it. It is negative
	* for large enough |eps| in the tails, where the expansion breaks
	* down.
	*
	* @param eps Size of the skew, from skew_eps.
	* @param w_par Velocity relative to the flow along B, in units of
	*        v_th.
	* @param w_sq Square of the velocity relative to the flow, in units of
	*        v_th^2.
	* @param shape Shape of the skew, from skew_shape.
	* @return The factor f / f_M.
	*/
	inline double skew_factor(const double eps, const double w_par,
		const double w_sq, const SkewShape& shape)
	{
		return 1.0 + eps * w_par
			* (shape.a + shape.b * w_sq + shape.c * w_sq * w_sq);
	}

	/**
	* Bound M on skew_factor used for rejection sampling (see
	* sample_bkg_velocity): M = 1 + |eps| shape.bound, which is at least the
	* factor for speeds up to 5 v_th. Faster draws (about 1 in 10^5 of a
	* Maxwellian) can exceed it; they are accepted with probability 1
	* instead of F / M, which is counted.
	*
	* @param eps Size of the skew.
	* @param shape Shape of the skew.
	* @return M.
	*/
	inline double skew_bound(const double eps, const SkewShape& shape)
	{
		return 1.0 + std::abs(eps) * shape.bound;
	}
	///@}
}
