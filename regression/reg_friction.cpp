/**
* @file reg_friction.cpp
*
* @brief Regression case: friction between impurities and a flowing plasma
*
* Tungsten ions start at rest in a uniform plasma that flows along the
* magnetic field. Collisions with the plasma drag them up to the flow
* speed. For a heavy impurity moving slowly compared with the background
* ions, the fluid (Spitzer) result is
*
*   <v_par>(t) = u (1 - exp(-t / tau_s)),
*
* with u the flow speed and tau_s the slowing-down time from collisions
* with the ions and the electrons. check_collisions.py reads the output,
* works out tau_s from the inputs saved in it, and compares.
*
* Everything about the case is set in this file: the plasma comes from
* the slab background (bkg_source = "slab"), which has nothing hard-coded
* about the plasma. Ionization and recombination are off, so the charge
* stays at its starting value.
*
* Usage: reg_friction <variant>, where <variant> is one of
*
*   z10        W10+, collision sub-cycling on (the defaults)
*   z10_nosub  W10+, sub-cycling off: shows what sub-cycling corrects
*   z5         W5+, to check the Z^2 scaling of the friction
*   z10_gpu    as z10, on the GPU
*
* The output is reg_friction_<variant>.nc in the working directory.
*/
#include <iostream>
#include <string>

#include "flan.h"
#include "flan_types.h"
#include "mpi.h"

/**
* Inputs of one variant of the case.
*
* @param variant Name of the variant, used in the case name.
* @param charge Charge state of the tungsten.
* @param tmax Length of the run (s), several slowing-down times.
* @param dt Time step (s).
* @param subcycle Whether the collisions are sub-cycled.
* @param gpu Whether to run on the GPU.
* @return The inputs.
*/
Inputs create_inputs(const std::string& variant, const int charge,
	const double tmax, const double dt, const bool subcycle, const bool gpu)
{
	Inputs inpts;

	inpts["case_name"] = "reg_friction_" + variant;
	inpts["use_gpu"] = gpu ? "cuda" : "off";

	// Plasma: uniform deuterium plasma flowing along B (+X). The flow is
	// under a fifth of the ion thermal speed (22 km/s), so the friction
	// is in its linear regime (to about 1%), and well above the
	// statistical noise of the mean impurity velocity.
	inpts["bkg_source"] = "slab";
	inpts["gkyl_ion_mass_amu"] = 2.014;
	inpts["gkyl_elec_mass_amu"] = 0.000548;
	inpts["slab_ne"] = 1e19;
	inpts["slab_te"] = 10.0;
	inpts["slab_ti"] = 10.0;
	inpts["slab_flow"] = 4000.0;
	inpts["slab_heat_flux_i"] = 0.0;
	inpts["slab_heat_flux_e"] = 0.0;
	inpts["slab_bmag"] = 1.0;
	inpts["slab_tmax"] = tmax;

	// Particles leaving one side of the box come back in the other
	inpts["min_xbound_type"] = "periodic";
	inpts["max_xbound_type"] = "periodic";
	inpts["min_ybound_type"] = "periodic";
	inpts["max_ybound_type"] = "periodic";
	inpts["min_zbound_type"] = "periodic";
	inpts["max_zbound_type"] = "periodic";

	// Tungsten, all starting at t = 0 in the middle of the box, at rest
	inpts["imp_num"] = 20000;
	inpts["imp_atom_num"] = 74;
	inpts["imp_mass_amu"] = 183.84;
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
	inpts["imp_temp_start_val"] = 0.0;

	// Time step and collisions. No thermal force: there is no heat flux.
	inpts["imp_time_step_opt"] = "constant";
	inpts["imp_time_step"] = dt;
	inpts["imp_collisions"] = "nanbu";
	inpts["imp_thermal_force"] = "off";
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

	// The run lengths are 5 to 6 slowing-down times. The time steps are
	// 1/24 of a gyration (W10+ gyrates in 1.2 us at 1 T, W5+ in 2.4 us).
	// They are long for the collisions on purpose: s0, Nanbu's s at the
	// thermal speed, is about 3 for the ions with W10+, where the model
	// without sub-cycling gives well under the fluid friction.
	if (variant == "z10")
		inpts = create_inputs(variant, 10, 40e-6, 5e-8, true, false);
	else if (variant == "z10_nosub")
		inpts = create_inputs(variant, 10, 40e-6, 5e-8, false, false);
	else if (variant == "z5")
		inpts = create_inputs(variant, 5, 150e-6, 1e-7, true, false);
	else if (variant == "z10_gpu")
		inpts = create_inputs(variant, 10, 40e-6, 5e-8, true, true);
	else
	{
		std::cerr << "Usage: reg_friction <variant>, with <variant> one of "
			"z10, z10_nosub, z5, z10_gpu\n";
		MPI_Finalize();
		return 1;
	}

	flan(inpts);

	MPI_Finalize();
	return 0;
}
