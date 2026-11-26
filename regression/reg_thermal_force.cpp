/**
* @file reg_thermal_force.cpp
*
* @brief Regression case: thermal force on impurities from the parallel
*        heat fluxes of the plasma
*
* Tungsten ions start at rest in a uniform plasma that doesn't flow but
* carries ion and electron heat fluxes along the magnetic field, as a
* plasma with a temperature gradient along the field does. Collisions
* then push the impurities against the heat flux, toward where the plasma
* would be hotter: the thermal force. They speed up until friction with
* the plasma cancels it, and from then on drift at
*
*   v_drift = -(F_i + F_e) / (m_z / tau_s),
*
* with the fluid forces F_i = beta_i Z^2 grad T_i and F_e = 0.71 Z^2
* grad T_e, beta_i Stangeby's mass-dependent coefficient (2.61 for W in
* D, 2.65 in the heavy limit, see below), grad T = q_par / kappa from Braginskii's conductivities, and
* tau_s the slowing-down time of reg_friction.cpp. Both the forces and the
* friction go as Z^2, so the drift doesn't depend on the charge.
* check_collisions.py reads the output, works all this out from the inputs
* saved in it, and compares.
*
* The heat fluxes are given as the size of the skew they put on the
* distributions, eps = q_par / (n T v_th). The ion one is 0.1, and the
* electron one is 0.0572, which makes the two temperature gradients
* equal (as for Ti = Te everywhere).
*
* Everything about the case is set in this file: the plasma comes from
* the slab background (bkg_source = "slab"). Ionization and recombination
* are off, so the charge stays at its starting value.
*
* Usage: reg_thermal_force <variant>, where <variant> is one of
*
*   z10        W10+, both heat fluxes, collision sub-cycling on
*   z10_ions   W10+, ion heat flux only: the 2.65 on its own
*   z10_elec   W10+, electron heat flux only: the 0.71 on its own.
*              Expected to miss the fluid drift by about 30%, see below.
*   z10_nosub  as z10 with sub-cycling off: shows what it corrects
*   z10_off    as z10 with imp_thermal_force = "off": no drift
*   z10_gpu    as z10, on the GPU
*   z5         W5+: the same drift as z10
*   c4_ions    C4+ starting at the ion temperature, ion heat flux only:
*              the thermal force on a light impurity (see below)
*
* The output is reg_thermal_force_<variant>.nc in the working directory.
*
* The cap on sub-steps is too small for the electrons here. With W10+ in
* this plasma and this time step, s0 (Nanbu's s at the thermal speed) is
* about 155 for collisions with electrons. The default cap of 128
* sub-steps leaves each at s0 = 1.2, far above the target of 0.03, and
* the electron thermal force comes out at about two thirds of its fluid
* value. Reaching the target would take about 5000 sub-steps per time
* step, or a z10_elec run of most of an hour. So:
*
* - z10_elec misses the fluid drift by about 30%. check_collisions.py
*   reports it as an expected miss and still holds it to what the
*   collision model should give with these settings.
* - z10 and z5 are within 10% of fluid all the same, since the electrons
*   carry only about a fifth of the thermal force, but z10 is close to
*   the edge.
*
* The light impurity. The heat-flux skew belongs to the background and is
* fitted so that a slow, heavy impurity gets 2.65 Z^2 grad T_i. A light
* impurity is not slow compared with the ions, and how the force changes
* with its mass comes only from the collisions themselves, so c4_ions
* compares Flan with Stangeby's mass-dependent coefficient,
*
*   beta_i = 3 (mu + 5 sqrt(2) Z^2 (1.1 mu^(5/2) - 0.35 mu^(3/2)) - 1)
*            / (2.6 - 2 mu + 5.4 mu^2),   mu = m_z / (m_z + m_i),
*
* about 2.05 Z^2 for C in D, against 2.65 Z^2 in the heavy limit. Its
* friction rate is that of two Maxwellians, and reg_friction c4 measures
* it in the same plasma.
*
* To check the 0.71 against fluid here, raise imp_coll_subcycle_cap for
* z10_elec to about 5000 and expect a long run. In the limit of small
* steps it is checked by the unit test ThermalForce.SmallStepLimit.
*/
#include <iostream>
#include <string>

#include "flan.h"
#include "flan_types.h"
#include "mpi.h"

/**
* The impurity of a variant.
*/
struct Impurity
{
	int atom_num {};		///< Atomic number
	double mass_amu {};		///< Mass (amu)
	double temp_start {};	///< Starting temperature (eV); 0 = at rest
	int num {};				///< Number of impurities followed
};

/// Tungsten, starting at rest
constexpr Impurity tungsten {74, 183.84, 0.0, 10000};

/// Carbon, starting at the ion temperature of the plasma (10 eV) so that
/// its velocity spread stays the same through the run. Its thermal speed
/// is 4 times tungsten's, so more of them are needed for the same
/// statistical error of the drift.
constexpr Impurity carbon {6, 12.011, 10.0, 20000};

/**
* Inputs of one variant of the case.
*
* @param variant Name of the variant, used in the case name.
* @param imp The impurity.
* @param charge Charge state of the impurity.
* @param tmax Length of the run (s), about 11 slowing-down times, so the
*        second half of the run is at force balance.
* @param dt Time step (s).
* @param eps_i Ion heat flux / (n Ti v_th,i).
* @param eps_e Electron heat flux / (n Te v_th,e).
* @param subcycle Whether the collisions are sub-cycled.
* @param thermal_force Whether the collisions include the thermal force.
* @return The inputs.
*/
Inputs create_inputs(const std::string& variant, const Impurity& imp,
	const int charge,
	const double tmax, const double dt, const double eps_i,
	const double eps_e, const bool subcycle, const bool thermal_force)
{
	Inputs inpts;

	inpts["case_name"] = "reg_thermal_force_" + variant;
	inpts["use_gpu"] = "off";

	// Plasma: uniform deuterium plasma at rest, with heat flowing along
	// B (+X)
	inpts["bkg_source"] = "slab";
	inpts["gkyl_ion_mass_amu"] = 2.014;
	inpts["gkyl_elec_mass_amu"] = 0.000548;
	inpts["slab_ne"] = 1e19;
	inpts["slab_te"] = 10.0;
	inpts["slab_ti"] = 10.0;
	inpts["slab_flow"] = 0.0;
	inpts["slab_heat_flux_i"] = eps_i;
	inpts["slab_heat_flux_e"] = eps_e;
	inpts["slab_bmag"] = 1.0;
	inpts["slab_tmax"] = tmax;

	// Particles leaving one side of the box come back in the other
	inpts["min_xbound_type"] = "periodic";
	inpts["max_xbound_type"] = "periodic";
	inpts["min_ybound_type"] = "periodic";
	inpts["max_ybound_type"] = "periodic";
	inpts["min_zbound_type"] = "periodic";
	inpts["max_zbound_type"] = "periodic";

	// Impurities, all starting at t = 0 in the middle of the box with no
	// mean velocity
	inpts["imp_num"] = imp.num;
	inpts["imp_atom_num"] = imp.atom_num;
	inpts["imp_mass_amu"] = imp.mass_amu;
	inpts["imp_init_charge"] = charge;
	inpts["imp_iz_recomb"] = "off";
	inpts["imp_tstart_opt"] = "single_value";
	inpts["imp_tstart_val"] = 0.0;
	inpts["imp_xstart_opt"] = "single_value";
	inpts["imp_xstart_val"] = 0.025;
	inpts["imp_ystart_opt"] = "single_value";
	inpts["imp_ystart_val"] = 0.0;
	inpts["imp_zstart_opt"] = "single_value";
	inpts["imp_zstart_val"] = 0.0;
	inpts["imp_temp_start_opt"] = "single_value";
	inpts["imp_temp_start_val"] = imp.temp_start;

	// Time step and collisions
	inpts["imp_time_step_opt"] = "constant";
	inpts["imp_time_step"] = dt;
	inpts["imp_collisions"] = "nanbu";
	inpts["imp_thermal_force"] = thermal_force ? "on" : "off";
	// The cap of 128 (the default) is too small for the electrons in
	// this plasma: see the note at the top of the file.
	inpts["imp_coll_subcycle_target"] = 0.03;
	inpts["imp_coll_subcycle_cap"] = subcycle ? 128 : 1;

	inpts["save_track"] = "off";
	inpts["print_interval"] = 10;

	return inpts;
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	const std::string variant {argc > 1 ? argv[1] : ""};
	Inputs inpts;

	// Skews of the two distributions. The electron one gives the same
	// temperature gradient as the ion one: eps_e / eps_i = (3.16 / 3.906)
	// / sqrt(2), from Braginskii's conductivities and collision times.
	constexpr double eps_i {0.1};
	constexpr double eps_e {0.0572};

	// Run lengths and time steps as in reg_friction.cpp, but twice as
	// long, to average the drift over the second half
	if (variant == "z10")
		inpts = create_inputs(variant, tungsten, 10, 80e-6, 5e-8,
			eps_i, eps_e, true, true);
	else if (variant == "z10_ions")
		inpts = create_inputs(variant, tungsten, 10, 80e-6, 5e-8,
			eps_i, 0.0, true, true);
	else if (variant == "z10_elec")
		inpts = create_inputs(variant, tungsten, 10, 80e-6, 5e-8,
			0.0, eps_e, true, true);
	else if (variant == "z10_nosub")
		inpts = create_inputs(variant, tungsten, 10, 80e-6, 5e-8,
			eps_i, eps_e, false, true);
	else if (variant == "z10_gpu")
	{
		inpts = create_inputs(variant, tungsten, 10, 80e-6, 5e-8,
			eps_i, eps_e, true, true);
		inpts["use_gpu"] = "cuda";
	}
	else if (variant == "z10_off")
		inpts = create_inputs(variant, tungsten, 10, 80e-6, 5e-8,
			eps_i, eps_e, true, false);
	else if (variant == "z5")
		inpts = create_inputs(variant, tungsten, 5, 300e-6, 1e-7,
			eps_i, eps_e, true, true);
	// C4+: 2.7 us slowing-down time, 0.2 us gyration. About 22
	// slowing-down times, so the drift is averaged over many of them.
	else if (variant == "c4_ions")
		inpts = create_inputs(variant, carbon, 4, 60e-6, 2e-8, eps_i, 0.0,
			true, true);
	else
	{
		std::cerr << "Usage: reg_thermal_force <variant>, with <variant> "
			"one of z10, z10_ions, z10_elec, z10_nosub, z10_off, z10_gpu, "
			"z5, c4_ions\n";
		MPI_Finalize();
		return 1;
	}

	flan(inpts);

	MPI_Finalize();
	return 0;
}
