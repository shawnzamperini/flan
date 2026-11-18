/**
* @file braginskii.h
* @brief Braginskii collision frequencies and parallel heat conduction, in
* the forms Hermes-3 uses (braginskii_collisions, braginskii_conduction).
*
* Used to calculate the parallel heat fluxes of backgrounds that don't
* store them (HERMES), so the thermal force can be included. Everything
* is in SI units, with temperatures in eV and densities in m^-3.
*/
#ifndef FLAN_BRAGINSKII_H
#define FLAN_BRAGINSKII_H

namespace Braginskii
{
	/// Constants (SI). qe and mp are the ones BOUT++ normalizes with (as
	/// in read_hermes.h); me and eps0 are CODATA 2018.
	inline constexpr double qe {1.602176634e-19};
	inline constexpr double mp {1.67262192369e-27};
	inline constexpr double me {9.1093837015e-31};
	inline constexpr double eps0 {8.8541878128e-12};

	/// Hermes-3's default parallel conduction coefficients: 3.9 for ions,
	/// and Braginskii's 3.16 divided by sqrt(2) for electrons, since
	/// Hermes-3's electron collision time is defined without the sqrt(2)
	/// that Braginskii's has.
	inline constexpr double kappa_coef_e_default {3.16 / 1.4142135623730951};
	inline constexpr double kappa_coef_i_default {3.9};

	double coulomb_log_ee(const double ne, const double te);

	double coulomb_log_ei(const double ne, const double ni, const double te,
		const double ti, const double z, const double aa);

	double coulomb_log_ii(const double n1, const double t1, const double z1,
		const double aa1, const double n2, const double t2, const double z2,
		const double aa2);

	double nu_ee(const double ne, const double te);

	double nu_ei(const double ne, const double ni, const double te,
		const double ti, const double z, const double aa);

	double nu_ii(const double n1, const double t1, const double z1,
		const double aa1, const double n2, const double t2, const double z2,
		const double aa2);

	double kappa_par(const double coef, const double n, const double t,
		const double nu, const double mass);

	double limit_kappa(const double kappa, const double grad_par_t,
		const double n, const double t, const double mass,
		const double alpha);
}

#endif
