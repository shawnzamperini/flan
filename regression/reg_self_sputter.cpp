/**
* @file reg_self_sputter.cpp
*
* @brief Regression cases: self-sputtering of tungsten on the slab wall
*
* Tungsten atoms are sputtered off a wall at y min of the slab background
* by the D+ flow onto it, as in reg_prompt_redep. With sput_self = "on",
* each W atom or ion landing on the wall sputters Y W atoms, and the
* particle carries on as one of them with its weight times Y
* (Sputtering::SelfSource). The weight carried into each generation g
* (W_g, with W_0 = S the D+ source) and the weight landing on the wall
* from it (L_g = f_g W_g) are tallied.
*
* With a constant yield Y0 the chain is simple: W_{g+1} = Y0 L_g, and if
* the atoms sputtered by W come back as often as those sputtered by D+
* (f_g = f), the self-sputtered rate is
*
*     T = sum_{g >= 1} W_g = S f Y0 / (1 - f Y0).
*
* Prompt redeposition doesn't depend on the atom's speed (p = Omega
* tau_iz), so the W-sputtered atoms, which are faster than the D-sputtered
* ones, should come back about as often. check_self_sputter.py checks
* these from the output.
*
* Usage: reg_self_sputter <variant>, where <variant> is one of
*
*   y05      Y0 = 0.5 (self_y05.h5): the weight halves each generation, so
*            from the fourth on it's below sput_self_min_weight (0.1 of
*            the source's) and Russian roulette is played.
*   y125     Y0 = 1.25 (self_y125.h5): the weight grows each generation,
*            but with f below 0.8 the chain still ends.
*   w_te<Te> The RustBCA W on W table (data/sputtering/W_on_W.h5) at
*            Te = Ti = <Te> eV, with full ionization and the sheath field
*            on ("brooks") at B 2 degrees to the wall. Printed by the
*            checker beside the yield at the energies the landing charges
*            fall through, not checked.
*
* The constant-yield runs are p1 of reg_prompt_redep (B = 2 T, n_e =
* 3e18 m^-3, Te = Ti = 20 eV, B at 0.01 degrees to the wall, no electric
* field, sheath, collisions or thermal force, ionizing no further than
* W+), so about 60% come back promptly and nearly all the rest are still
* gyrating out in the plasma when they're stopped.
*
* The tables self_y05.h5 and self_y125.h5 are made by
* make_self_tables.py in the working directory. The D on W and W on W
* tables are read from data/sputtering/ in Flan's directory
* (FLAN_DATA_DIR, set by CMakeLists.txt), and the OpenADAS rates from the
* directory in the environment variable OPENADAS_ROOT. The output is
* reg_self_sputter_<variant>.nc in the working directory.
*/
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

#include "flan.h"
#include "flan_types.h"
#include "mpi.h"

#ifndef FLAN_DATA_DIR
#define FLAN_DATA_DIR "../data"
#endif

/**
* Inputs shared by every variant: a straight slab with a sputtered wall at
* y min, B = 2 T along +X turned 0.01 degrees toward the wall, D plasma
* flowing along B into the wall at 20 eV, and 20000 sputtered W atoms
* followed without collisions, ionizing no further than W+, with
* self-sputtering on.
*
* @param variant Name of the variant, used in the case name.
* @return The inputs.
*/
Inputs base_inputs(const std::string& variant)
{
	Inputs inpts;

	inpts["case_name"] = "reg_self_sputter_" + variant;
	inpts["use_gpu"] = "off";

	// Plasma
	inpts["bkg_source"] = "slab";
	inpts["gkyl_ion_mass_amu"] = 2.014;
	inpts["gkyl_elec_mass_amu"] = 0.000548;
	inpts["slab_ne"] = 3e18;
	inpts["slab_te"] = 20.0;
	inpts["slab_ti"] = 20.0;
	inpts["slab_flow"] = 2e4;
	inpts["slab_heat_flux_i"] = 0.0;
	inpts["slab_heat_flux_e"] = 0.0;
	inpts["slab_bmag"] = 2.0;
	inpts["slab_efield_y"] = 0.0;
	inpts["slab_efield_y_rate"] = 0.0;
	inpts["slab_bgrad_y"] = 0.0;
	inpts["slab_major_radius"] = 0.0;

	// The wall, with B nearly parallel to it
	inpts["slab_wall"] = "ymin";
	inpts["slab_b_angle"] = 0.01;

	// Periodic along the wall, absorbing at the far side
	inpts["min_xbound_type"] = "periodic";
	inpts["max_xbound_type"] = "periodic";
	inpts["min_ybound_type"] = "absorbing";
	inpts["max_ybound_type"] = "absorbing";
	inpts["min_zbound_type"] = "periodic";
	inpts["max_zbound_type"] = "periodic";
	inpts["tbound_type"] = "periodic";

	// Sputtered tungsten
	inpts["imp_num"] = 20000;
	inpts["imp_atom_num"] = 74;
	inpts["imp_mass_amu"] = 183.84;
	inpts["imp_source_opt"] = "sputter";
	inpts["sput_table_file"] = std::string(FLAN_DATA_DIR)
		+ "/sputtering/D_on_W.h5";
	inpts["sput_angle_val"] = 60.0;

	// Self-sputtering
	inpts["sput_self"] = "on";
	inpts["sput_self_min_weight"] = 0.1;
	inpts["sput_self_max_gen"] = 1000;

	// Atomic physics, ionizing no further than W+
	inpts["imp_iz_recomb"] = "on";
	inpts["imp_iz_max_charge"] = 1;
	const char* adas {std::getenv("OPENADAS_ROOT")};
	inpts["openadas_root"] = std::string(adas != nullptr ? adas : "undefined");
	inpts["openadas_year"] = 50;

	inpts["imp_time_step_opt"] = "constant";
	inpts["imp_time_step"] = 5e-9;
	inpts["imp_collisions"] = "off";
	inpts["imp_thermal_force"] = "off";
	inpts["sheath_model"] = "off";

	inpts["save_track"] = "off";
	inpts["print_interval"] = 10;

	return inpts;
}

/**
* Sets how long each atom is followed (from its own sputtering) and the
* background's length: 7 mean ionization times of W0 plus 2 gyrations of
* W+, using a rough S0 of 3.1e-13 m^3/s at 20 eV, as in reg_prompt_redep.
*
* @param inpts Inputs to change.
* @param ne Electron density (m^-3).
* @param B Field strength (T).
*/
void set_times(Inputs& inpts, const double ne, const double B)
{
	constexpr double S0 {3.1e-13};
	constexpr double m {183.84 * 1.661e-27};
	constexpr double e {1.602e-19};
	const double t_iz {1.0 / (ne * S0)};
	const double t_gyro {2.0 * 3.141592653589793 * m / (e * B)};
	const double t_max {7.0 * t_iz + 2.0 * t_gyro};
	inpts["slab_ne"] = ne;
	inpts["slab_bmag"] = B;
	inpts["imp_max_time"] = t_max;
	inpts["slab_tmax"] = t_max;
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	if (std::getenv("OPENADAS_ROOT") == nullptr)
	{
		std::cerr << "Set OPENADAS_ROOT to the directory holding the "
			"OpenADAS files (scd50/, acd50/).\n";
		MPI_Finalize();
		return 1;
	}

	const std::string variant {argc > 1 ? argv[1] : ""};
	Inputs inpts {base_inputs(variant)};

	// Constant yields, from make_self_tables.py
	if (variant == "y05" || variant == "y125")
	{
		const std::string table {"self_" + variant + ".h5"};
		if (!std::filesystem::exists(table))
		{
			std::cerr << table << " isn't in the working directory. Make it "
				"with make_self_tables.py first.\n";
			MPI_Finalize();
			return 1;
		}
		set_times(inpts, 3e18, 2.0);
		inpts["sput_self_table_file"] = table;
	}

	// RustBCA's W on W, at Te = Ti = <Te> eV with full ionization and the
	// sheath field on, B at 2 degrees to the wall
	else if (variant.starts_with("w_te"))
	{
		const double te {std::stod(variant.substr(4))};
		set_times(inpts, 3e18, 2.0);
		inpts["slab_te"] = te;
		inpts["slab_ti"] = te;
		inpts["slab_flow"] = std::sqrt(2.0 * te * 1.602e-19
			/ (2.014 * 1.661e-27));
		inpts["slab_b_angle"] = 2.0;
		inpts["imp_iz_max_charge"] = -1;
		inpts["sheath_model"] = "brooks";
		inpts["sput_self_table_file"] = std::string(FLAN_DATA_DIR)
			+ "/sputtering/W_on_W.h5";
	}
	else
	{
		std::cerr << "Usage: reg_self_sputter <variant>, with <variant> one "
			"of y05, y125, w_te<Te>\n";
		MPI_Finalize();
		return 1;
	}

	flan(inpts);

	MPI_Finalize();
	return 0;
}
