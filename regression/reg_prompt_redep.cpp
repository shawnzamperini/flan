/**
* @file reg_prompt_redep.cpp
*
* @brief Regression cases: prompt redeposition of sputtered tungsten
*        without a sheath (Fussmann's model)
*
* Tungsten atoms are sputtered off a wall at y min of the slab background
* with a cosine distribution of directions, ionize, and gyrate. With the
* magnetic field parallel to the wall and no electric field, an atom that
* ionizes after a path of x gyroradii (of the singly charged ion at the
* atom's speed) comes back to the wall in its first gyration with
* probability 1 / (1 + x^2). Averaged over the exponential distribution of
* ionization lengths, the fraction that comes back promptly depends only
* on p = e B / (m n_e S), the mean ionization length over the gyroradius
* (S the ionization rate coefficient of the atom):
*
*     <r>(p) = (1 / p) int_0^inf exp(-x / p) / (1 + x^2) dx.
*
* (See regression/fussmann.py.) check_redeposition.py compares Flan's
* prompt fraction (returns within one gyroperiod of the singly charged ion
* after first ionizing, Sputtering::Deposition::prompt_total, over the
* source) with <r>(p) for each run.
*
* Usage: reg_prompt_redep <variant>, where <variant> is one of
*
*   p01, p03, p1, p3, p10  B = 2 T and n_e from 3e19 to 3e17 m^-3, so p
*                          from about 0.1 to 10 (the names are rough; the
*                          checker works out p from the inputs)
*   b4                     B = 4 T and n_e = 6e18 m^-3: p about 1 again,
*                          checking that it's B / n_e that matters
*   full                   As p1, with ionization past W+ allowed. Not a
*                          check of the model (which leaves further
*                          ionization out); the checker prints it to show
*                          how much it matters.
*   drop0                  As p1 with the sheath field on ("brooks") but a
*                          drop of 1e-4 Te (sheath_drop): Fussmann's case
*                          again, checked the same way.
*   brooks_a<angle>,       As p1 with the sheath field on and B at <angle>
*   borodkina_a<angle>     degrees to the wall (2, 5, 10, 20 in the
*                          README): how much the sheath raises the prompt
*                          fraction, and how much that depends on the
*                          profile. Printed, not checked.
*   <variant>_gpu          Any of the above without the sheath (p01 to
*                          p10, b4, full), on the GPU: the source, the
*                          ionization with its charge limit and prompt
*                          window, and the deposition tallies run there,
*                          and are checked the same way.
*   brooks_ref             Brooks' reference case (1990, Table 1): D-T at
*                          75 eV, 1e20 m^-3, B = 5 T at 87 degrees from the
*                          normal, his f_D = 0.25 and 3 Te drop
*                          (sheath_fd, sheath_drop), collisions on, full
*                          ionization. The checker compares the
*                          redeposited W's charge, impact energy and
*                          angle, and time to return with his.
*
* The model's assumptions, and how the runs meet them:
*
* - Field parallel to the wall. The plasma has to reach the wall to
*   sputter it, so the field is at 0.01 degrees to it (slab_b_angle). In
*   one gyration an ion then moves about 1e-3 of a gyroradius toward or
*   away from the wall along B, which is negligible.
* - No electric field, collisions or thermal force: all off, the sheath
*   field included (sheath_model = "off").
* - No further ionization in the first gyration. W+ would ionize again
*   about 2 pi (S_1 / S_0) / p times per gyration, close to once at p = 1,
*   so the runs stop ionization at W+ (imp_iz_max_charge = 1).
* - The atom's speed doesn't matter, since p doesn't depend on it, so the
*   normal Thompson energy distribution is used.
*
* Each run follows an atom for 7 mean ionization times plus 2 gyrations,
* long enough for nearly every atom to ionize and finish its first
* gyration. An atom that hasn't by then would have ionized so far out
* that it comes back with probability below 1e-3.
*
* Everything about a case is set in this file through the slab_* options;
* the plasma is uniform at Te = Ti = 20 eV. Two files come from outside:
* the D on W yield table, data/sputtering/D_on_W.h5 in Flan's directory
* (FLAN_DATA_DIR, set by CMakeLists.txt), and the OpenADAS rates, from the
* directory in the environment variable OPENADAS_ROOT. The output is
* reg_prompt_redep_<variant>.nc in the working directory.
*/
#include <cmath>
#include <cstdlib>
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
* flowing along B into the wall, and 20000 sputtered W atoms followed
* without collisions, ionizing no further than W+.
*
* @param variant Name of the variant, used in the case name.
* @return The inputs.
*/
Inputs base_inputs(const std::string& variant)
{
	Inputs inpts;

	inpts["case_name"] = "reg_prompt_redep_" + variant;
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
	// The D+ impacts from the formula at sput_angle_val, not the sheath
	// model (the default), which some variants have off: this case
	// checks what the atoms do after they're sputtered
	inpts["sput_angle_val"] = 60.0;
	inpts["sput_impact_opt"] = "formula";

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
* Sets how long atoms are followed (and the background's length): 7 mean
* ionization times of W0 plus 2 gyrations of W+, using a rough S0 of
* 3.1e-13 m^3/s at 20 eV. The checker doesn't depend on this value.
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

	const std::string name {argc > 1 ? argv[1] : ""};
	Inputs inpts {base_inputs(name)};

	// <variant>_gpu: the same case on the GPU
	std::string variant {name};
	if (variant.ends_with("_gpu"))
	{
		variant.resize(variant.size() - 4);
		inpts["use_gpu"] = "cuda";
	}

	if (variant == "p01") set_times(inpts, 3e19, 2.0);
	else if (variant == "p03") set_times(inpts, 1e19, 2.0);
	else if (variant == "p1") set_times(inpts, 3e18, 2.0);
	else if (variant == "p3") set_times(inpts, 1e18, 2.0);
	else if (variant == "p10") set_times(inpts, 3e17, 2.0);
	else if (variant == "b4") set_times(inpts, 6e18, 4.0);
	else if (variant == "full")
	{
		set_times(inpts, 3e18, 2.0);
		inpts["imp_iz_max_charge"] = -1;
	}

	// With the sheath: a drop so small that it's Fussmann's case again
	else if (variant == "drop0")
	{
		set_times(inpts, 3e18, 2.0);
		inpts["sheath_model"] = "brooks";
		inpts["sheath_drop"] = 1e-4;
	}

	// With the sheath, as p1 with B at 2 to 20 degrees to the wall, for
	// either profile: <model>_a<angle>
	else if (variant.starts_with("brooks_a")
		|| variant.starts_with("borodkina_a"))
	{
		const auto u {variant.find("_a")};
		set_times(inpts, 3e18, 2.0);
		inpts["sheath_model"] = variant.substr(0, u);
		inpts["slab_b_angle"] = std::stod(variant.substr(u + 2));
	}

	// Brooks' reference case (Phys. Fluids B 2, 1858, 1990, Table 1): D-T
	// at Te = Ti = 75 eV and 1e20 m^-3, B = 5 T at 87 degrees from the
	// normal, his f_D = 0.25 and drop of 3 Te, collisions on, and
	// ionization as far as it goes
	else if (variant == "brooks_ref")
	{
		inpts["gkyl_ion_mass_amu"] = 2.5;
		inpts["slab_ne"] = 1e20;
		inpts["slab_te"] = 75.0;
		inpts["slab_ti"] = 75.0;
		inpts["slab_bmag"] = 5.0;
		inpts["slab_flow"] = std::sqrt(150.0 * 1.602e-19 / (2.5 * 1.661e-27));
		inpts["slab_b_angle"] = 3.0;
		inpts["sheath_model"] = "brooks";
		inpts["sheath_fd"] = 0.25;
		inpts["sheath_drop"] = 3.0;
		inpts["imp_iz_max_charge"] = -1;
		inpts["imp_collisions"] = "nanbu";
		inpts["imp_max_time"] = 1e-5;
		inpts["slab_tmax"] = 1e-5;
		inpts["imp_time_step"] = 1e-9;
	}
	else
	{
		std::cerr << "Usage: reg_prompt_redep <variant>, with <variant> one "
			"of p01, p03, p1, p3, p10, b4, full, drop0, brooks_a<angle>, "
			"borodkina_a<angle>, brooks_ref, or one without the sheath "
			"followed by _gpu\n";
		MPI_Finalize();
		return 1;
	}

	// The sheath field isn't on the GPU yet
	if (variant != name
		&& std::get<std::string>(inpts["sheath_model"]) != "off")
	{
		std::cerr << "The _gpu variants are for the cases without the "
			"sheath.\n";
		MPI_Finalize();
		return 1;
	}

	flan(inpts);

	MPI_Finalize();
	return 0;
}
