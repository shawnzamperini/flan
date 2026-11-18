/**
* @file braginskii.cpp
* @brief Braginskii collision frequencies and parallel heat conduction.
*
* These follow the forms in Hermes-3's braginskii_collisions and
* braginskii_conduction components (written here from the formulas, not
* the code), so that a heat flux calculated from a HERMES run's n and T is
* the one the run itself had. The Coulomb logarithms are the NRL
* Formulary's with densities in m^-3. Temperatures are in eV and
* densities in m^-3 throughout.
*
* Hermes-3 floors with a smooth softFloor; plain max() is used here,
* which only differs right at the floor.
*/

#include <algorithm>
#include <cmath>
#include <numbers>

#include "braginskii.h"

namespace Braginskii
{
	namespace
	{
		/// Temperature (eV) and density (m^-3) floors of the Coulomb
		/// logarithms and thermal speeds
		constexpr double t_floor {0.1};
		constexpr double n_floor {1e10};

		constexpr double pi {std::numbers::pi};

		/// Collision frequencies are proportional to this
		/// (qe^4 / (eps0 m)^2), for a particle of mass m
		double prefactor(const double m)
		{
			const double em {eps0 * m};
			return qe * qe * qe * qe / (em * em);
		}
	}

	/**
	* Electron-electron Coulomb logarithm (NRL Formulary):
	* 30.4 - ln(ne) / 2 + 5/4 ln(Te) - sqrt(1e-5 + (ln(Te) - 2)^2 / 16),
	* with Te floored at 0.1 eV and ne at 1e10 m^-3.
	*
	* @param ne Electron density (m^-3).
	* @param te Electron temperature (eV).
	* @return ln Lambda_ee (not floored; nu_ee floors it at 1).
	*/
	double coulomb_log_ee(const double ne, const double te)
	{
		const double t {std::max(te, t_floor)};
		const double n {std::max(ne, n_floor)};
		const double lt {std::log(t)};
		return 30.4 - 0.5 * std::log(n) + 1.25 * lt
			- std::sqrt(1e-5 + (lt - 2.0) * (lt - 2.0) / 16.0);
	}

	/**
	* Electron-ion Coulomb logarithm (NRL Formulary), in four ranges of Te:
	*   - 10 if Te < 0.1 eV or either density is under 1e10 m^-3
	*   - 23 - ln(ni) / 2 + 3/2 ln(Ti) - ln(Z^2 A) if Te < Ti me / mi
	*   - 30 - ln(ne) / 2 - ln(Z) + 3/2 ln(Te) if Te < e^2 Z^2 eV
	*   - 31 - ln(ne) / 2 + ln(Te) otherwise.
	* The third range ends at e^2 Z^2 (about 7.4 Z^2 eV), as in Hermes-3,
	* where the two expressions meet, rather than the NRL's 10 Z^2.
	*
	* @param ne Electron density (m^-3).
	* @param ni Ion density (m^-3).
	* @param te Electron temperature (eV).
	* @param ti Ion temperature (eV).
	* @param z Ion charge.
	* @param aa Ion mass (amu, in proton masses).
	* @return ln Lambda_ei (not floored; nu_ei floors it at 1).
	*/
	double coulomb_log_ei(const double ne, const double ni, const double te,
		const double ti, const double z, const double aa)
	{
		const double me_mi {me / (mp * aa)};
		if (te < t_floor || ni < n_floor || ne < n_floor) return 10.0;
		if (te < ti * me_mi)
		{
			return 23.0 - 0.5 * std::log(ni) + 1.5 * std::log(ti)
				- std::log(z * z * aa);
		}
		if (te < std::exp(2.0) * z * z)
		{
			return 30.0 - 0.5 * std::log(ne) - std::log(z)
				+ 1.5 * std::log(te);
		}
		return 31.0 - 0.5 * std::log(ne) + std::log(te);
	}

	/**
	* Ion-ion Coulomb logarithm (NRL Formulary, mixed ion-ion):
	* 29.91 - ln(Z1 Z2 (A1 + A2) / (A1 T2 + A2 T1)
	* sqrt(n1 Z1^2 / T1 + n2 Z2^2 / T2)), with temperatures floored at
	* 0.1 eV and densities at 1e10 m^-3. For like ions, pass the same
	* species twice.
	*
	* @param n1 Density of species 1 (m^-3).
	* @param t1 Temperature of species 1 (eV).
	* @param z1 Charge of species 1.
	* @param aa1 Mass of species 1 (amu).
	* @param n2 Density of species 2 (m^-3).
	* @param t2 Temperature of species 2 (eV).
	* @param z2 Charge of species 2.
	* @param aa2 Mass of species 2 (amu).
	* @return ln Lambda_12 (not floored; nu_ii floors it at 1).
	*/
	double coulomb_log_ii(const double n1, const double t1, const double z1,
		const double aa1, const double n2, const double t2, const double z2,
		const double aa2)
	{
		const double ta {std::max(t1, t_floor)};
		const double tb {std::max(t2, t_floor)};
		const double na {std::max(n1, n_floor)};
		const double nb {std::max(n2, n_floor)};
		return 29.91 - std::log(z1 * z2 * (aa1 + aa2) / (aa1 * tb + aa2 * ta)
			* std::sqrt(na * z1 * z1 / ta + nb * z2 * z2 / tb));
	}

	/**
	* Electron-electron collision frequency:
	* nu_ee = 2 qe^4 ne lnL / (3 (2 pi v^2)^(3/2) (eps0 me)^2), with
	* v^2 = 2 e Te / me, Te floored at 0.1 eV and lnL at 1.
	*
	* @param ne Electron density (m^-3).
	* @param te Electron temperature (eV).
	* @return nu_ee (s^-1).
	*/
	double nu_ee(const double ne, const double te)
	{
		const double v2 {2.0 * std::max(te, t_floor) * qe / me};
		const double lnl {std::max(coulomb_log_ee(ne, te), 1.0)};
		return 2.0 * prefactor(me) * std::max(ne, 0.0) * lnl
			/ (3.0 * std::pow(2.0 * pi * v2, 1.5));
	}

	/**
	* Electron-ion collision frequency (electrons colliding with ions):
	* nu_ei = qe^4 Z^2 ni lnL (1 + me / mi)
	* / (3 (pi (ve^2 + vi^2))^(3/2) (eps0 me)^2), with v^2 = 2 e T / m,
	* temperatures floored at 0.1 eV and lnL at 1.
	*
	* @param ne Electron density (m^-3).
	* @param ni Ion density (m^-3).
	* @param te Electron temperature (eV).
	* @param ti Ion temperature (eV).
	* @param z Ion charge.
	* @param aa Ion mass (amu).
	* @return nu_ei (s^-1).
	*/
	double nu_ei(const double ne, const double ni, const double te,
		const double ti, const double z, const double aa)
	{
		const double mi {mp * aa};
		const double ve2 {2.0 * std::max(te, t_floor) * qe / me};
		const double vi2 {2.0 * std::max(ti, t_floor) * qe / mi};
		const double lnl {std::max(coulomb_log_ei(ne, ni, te, ti, z, aa),
			1.0)};
		return prefactor(me) * z * z * std::max(ni, 0.0) * lnl
			* (1.0 + me / mi) / (3.0 * std::pow(pi * (ve2 + vi2), 1.5));
	}

	/**
	* Collision frequency of ion species 1 with ion species 2:
	* nu_12 = (Z1 Z2)^2 qe^4 n2 lnL (1 + m1 / m2)
	* / (3 (pi (v1^2 + v2^2))^(3/2) (eps0 m1)^2), with v^2 = 2 e T / m,
	* temperatures floored at 0.1 eV, n2 at 1e10 m^-3 and lnL at 1. For
	* like-ion collisions, pass the same species twice.
	*
	* @param n1 Density of species 1 (m^-3).
	* @param t1 Temperature of species 1 (eV).
	* @param z1 Charge of species 1.
	* @param aa1 Mass of species 1 (amu).
	* @param n2 Density of species 2 (m^-3).
	* @param t2 Temperature of species 2 (eV).
	* @param z2 Charge of species 2.
	* @param aa2 Mass of species 2 (amu).
	* @return nu_12 (s^-1).
	*/
	double nu_ii(const double n1, const double t1, const double z1,
		const double aa1, const double n2, const double t2, const double z2,
		const double aa2)
	{
		const double m1 {mp * aa1};
		const double m2 {mp * aa2};
		const double v1 {2.0 * std::max(t1, t_floor) * qe / m1};
		const double v2 {2.0 * std::max(t2, t_floor) * qe / m2};
		const double lnl {std::max(coulomb_log_ii(n1, t1, z1, aa1, n2, t2,
			z2, aa2), 1.0)};
		const double zz {z1 * z2};
		return zz * zz * prefactor(m1) * std::max(n2, n_floor) * lnl
			* (1.0 + m1 / m2) / (3.0 * std::pow(pi * (v1 + v2), 1.5));
	}

	/**
	* Parallel heat conductivity, kappa = coef n (e T) tau / m times e, so
	* that q_par = -kappa grad_par T with T in eV gives W/m^2 (one e makes
	* the pressure J/m^3, the other the gradient J/m). This is Hermes-3's
	* kappa = coef P tau / AA in SI units: tau = 1 / nu is the species'
	* collision time from all of its collisions, and the pressure is
	* floored at 0. With coef = 3.16 and only electron-ion collisions it's
	* Spitzer's kappa_e, about 3.1e4 Te^2.5 / lnL W / m / eV.
	*
	* @param coef Conduction coefficient (e.g., kappa_coef_e_default).
	* @param n Density of the species (m^-3).
	* @param t Temperature of the species (eV).
	* @param nu Sum of the species' collision frequencies (s^-1).
	* @param mass Mass of the species (kg).
	* @return kappa_par (W / m / eV), 0 if there's no pressure or no
	*         collisions.
	*/
	double kappa_par(const double coef, const double n, const double t,
		const double nu, const double mass)
	{
		const double p {n * t};
		if (!(p > 0.0) || !(nu > 0.0)) return 0.0;
		return coef * p * qe * qe / (mass * nu);
	}

	/**
	* Applies Hermes-3's optional flux limiter to a conductivity:
	* kappa / (1 + |q_SH / q_fl|), with q_SH = kappa grad_par T and the
	* free-streaming limit q_fl = alpha n e T sqrt(e T / m). This keeps the
	* heat flux from exceeding alpha times the free-streaming flux where
	* the plasma is hot and collisionless.
	*
	* @param kappa Conductivity without the limit (W / m / eV).
	* @param grad_par_t Parallel gradient of T (eV / m).
	* @param n Density of the species (m^-3).
	* @param t Temperature of the species (eV).
	* @param mass Mass of the species (kg).
	* @param alpha Flux limiter. Zero or negative means no limit.
	* @return The limited conductivity.
	*/
	double limit_kappa(const double kappa, const double grad_par_t,
		const double n, const double t, const double mass,
		const double alpha)
	{
		if (!(alpha > 0.0)) return kappa;
		const double q_fl {alpha * std::max(n, 0.0) * qe * std::max(t, 0.0)
			* std::sqrt(qe * std::max(t, 0.0) / mass)};
		const double q_sh {kappa * grad_par_t};
		if (!(q_fl > 0.0)) return q_sh == 0.0 ? kappa : 0.0;
		return kappa / (1.0 + std::abs(q_sh / q_fl));
	}
}
