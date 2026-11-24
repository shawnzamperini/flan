/**
* @file reg_sheath_impacts.cpp
*
* @brief Regression cases: the sputtered source from the D+ impact
*        energies and angles of the sheath model (sput_impact_opt)
*
* D+ flows onto a wall at y min of the slab background, with the field at
* an angle theta to the wall's normal. With sput_impact_opt = "sheath",
* Flan pushes D+ test ions through the sheath field (Sheath::ImpactTable)
* and averages the yield over the energies and angles they hit the wall
* with; with "formula" every ion hits with the energy the sheath drop and
* the temperatures give, at sput_angle_val. The source is all that's
* checked, so only 100 W atoms are followed, without atomic physics, for
* a few time steps.
*
* check_sheath_impacts.py checks, from the output:
*
*  - With a constant yield, the two options give the same source.
*  - With a yield linear in energy and flat in angle, the yield averaged
*    over the impacts is the slope times their mean energy, and that is
*    the drop plus the mean energy the ions enter the sheath with
*    (weighted by their share of the flux), since the field does work
*    only through the drop.
*  - At theta = 0 the field is along B, so only the speed along the
*    normal changes, and the angles they hit at follow from the starting
*    Maxwellian and the drop.
*  - In Borodkina and Tsvetkov's cases (reg_sheath_angles), the mean
*    angle matches the one-dimensional sampler at those exact plasma
*    values, which shows the error of the table's grid.
*
* and prints the source from RustBCA's D on W table with both options at
* a few Te and field angles.
*
* Usage: reg_sheath_impacts <variant>. A variant is a list of parts
* separated by "_": first the table and the plasma,
*
*   const    Y = 0.02 (impact_const.h5), Brooks' sheath, n = 1e19 m^-3,
*            B = 2 T
*   lin      Y = 1e-3 E / eV, flat in angle (impact_lin.h5), as const
*   bor      Y = 0.02, Borodkina's sheath and her cases (with a45 or a60:
*            Te = Ti = 30 eV, n = 1e18 m^-3 and B = 5 T, or n = 1e20 m^-3
*            and B = 3 T)
*   dw       RustBCA's D on W (data/sputtering/D_on_W.h5), as const
*
* then any of
*
*   a<deg>   theta, the field's angle to the normal (degrees, default 60)
*   te<Te>   Te in eV (default 20)
*   ti<r>    Ti / Te (default 1), moved to the nearest point of the
*            impact table's grid in ln(Ti / Te), so the table is exact
*   formula  sput_impact_opt = "formula" (default "sheath")
*
* e.g. lin_a85_ti2 or dw_te50_a85_formula. The tables impact_const.h5 and
* impact_lin.h5 are made by make_impact_tables.py in the working
* directory; D_on_W.h5 is read from data/sputtering/ in Flan's directory
* (FLAN_DATA_DIR, set by CMakeLists.txt). The output is
* reg_sheath_impacts_<variant>.nc in the working directory.
*/
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "flan.h"
#include "flan_types.h"
#include "mpi.h"

#ifndef FLAN_DATA_DIR
#define FLAN_DATA_DIR "../data"
#endif

/**
* Inputs shared by every variant: a straight slab with a sputtered wall at
* y min, D plasma flowing along B into it at the sound speed, and 100 W
* atoms followed for 20 steps without atomic physics or collisions.
*
* @param variant Name of the variant, used in the case name.
* @return The inputs.
*/
Inputs base_inputs(const std::string& variant)
{
	Inputs inpts;

	inpts["case_name"] = "reg_sheath_impacts_" + variant;
	inpts["use_gpu"] = "off";

	// Plasma
	inpts["bkg_source"] = "slab";
	inpts["gkyl_ion_mass_amu"] = 2.014;
	inpts["gkyl_elec_mass_amu"] = 0.000548;
	inpts["slab_ne"] = 1e19;
	inpts["slab_te"] = 20.0;
	inpts["slab_ti"] = 20.0;
	inpts["slab_heat_flux_i"] = 0.0;
	inpts["slab_heat_flux_e"] = 0.0;
	inpts["slab_bmag"] = 2.0;
	inpts["slab_efield_y"] = 0.0;
	inpts["slab_efield_y_rate"] = 0.0;
	inpts["slab_bgrad_y"] = 0.0;
	inpts["slab_major_radius"] = 0.0;
	inpts["slab_wall"] = "ymin";

	inpts["min_xbound_type"] = "periodic";
	inpts["max_xbound_type"] = "periodic";
	inpts["min_ybound_type"] = "absorbing";
	inpts["max_ybound_type"] = "absorbing";
	inpts["min_zbound_type"] = "periodic";
	inpts["max_zbound_type"] = "periodic";
	inpts["tbound_type"] = "periodic";

	// Sputtered tungsten, only a few followed: the source is what's
	// checked
	inpts["imp_num"] = 100;
	inpts["imp_atom_num"] = 74;
	inpts["imp_mass_amu"] = 183.84;
	inpts["imp_source_opt"] = "sputter";
	inpts["sput_angle_val"] = 60.0;
	inpts["sheath_model"] = "brooks";
	inpts["sput_impact_opt"] = "sheath";

	inpts["imp_iz_recomb"] = "off";
	inpts["imp_time_step_opt"] = "constant";
	inpts["imp_time_step"] = 5e-9;
	inpts["imp_max_time"] = 1e-7;
	inpts["slab_tmax"] = 1e-7;
	inpts["imp_collisions"] = "off";
	inpts["imp_thermal_force"] = "off";

	inpts["save_track"] = "off";
	inpts["print_interval"] = 10;

	return inpts;
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	const std::string variant {argc > 1 ? argv[1] : ""};
	Inputs inpts {base_inputs(variant)};

	std::vector<std::string> parts {};
	{
		std::istringstream ss {variant};
		std::string part {};
		while (std::getline(ss, part, '_')) parts.push_back(part);
	}

	auto usage = [&]()
	{
		std::cerr << "Usage: reg_sheath_impacts <variant>, with <variant> "
			"const, lin, bor or dw followed by any of _a<deg>, _te<Te>, "
			"_ti<Ti/Te>, _formula (see reg_sheath_impacts.cpp)\n";
		MPI_Finalize();
		return 1;
	};
	if (parts.empty()) return usage();

	// The table and the plasma
	const std::string kind {parts[0]};
	double theta {60.0};
	double te {20.0};
	double ti_te {1.0};
	if (kind == "const" || kind == "bor")
		inpts["sput_table_file"] = std::string {"impact_const.h5"};
	else if (kind == "lin")
		inpts["sput_table_file"] = std::string {"impact_lin.h5"};
	else if (kind == "dw")
	{
		inpts["sput_table_file"] = std::string(FLAN_DATA_DIR)
			+ "/sputtering/D_on_W.h5";
	}
	else return usage();

	try
	{
		for (std::size_t i {1}; i < parts.size(); ++i)
		{
			const std::string& p {parts[i]};
			if (p == "formula") inpts["sput_impact_opt"] = "formula";
			else if (p.starts_with("a")) theta = std::stod(p.substr(1));
			else if (p.starts_with("te")) te = std::stod(p.substr(2));
			else if (p.starts_with("ti")) ti_te = std::stod(p.substr(2));
			else return usage();
		}
	}
	catch (const std::exception&) {return usage();}

	const auto table {std::get<std::string>(inpts["sput_table_file"])};
	if (!std::filesystem::exists(table))
	{
		std::cerr << table << " isn't there. Make impact_const.h5 and "
			"impact_lin.h5 with make_impact_tables.py first.\n";
		MPI_Finalize();
		return 1;
	}

	// Borodkina and Tsvetkov's cases, as in reg_sheath_angles
	if (kind == "bor")
	{
		inpts["sheath_model"] = "borodkina";
		te = 30.0;
		ti_te = 1.0;
		if (theta == 45.0)
		{
			inpts["slab_ne"] = 1e18;
			inpts["slab_bmag"] = 5.0;
		}
		else if (theta == 60.0)
		{
			inpts["slab_ne"] = 1e20;
			inpts["slab_bmag"] = 3.0;
		}
		else return usage();
	}

	// Ti / Te on the impact table's grid, so its samples are for this
	// plasma exactly (energies scale with Te, and theta is checked on the
	// grid's 0.5 degree points)
	const double step {0.1};
	ti_te = std::exp(step * std::round(std::log(ti_te) / step));
	inpts["slab_te"] = te;
	inpts["slab_ti"] = te * ti_te;
	inpts["slab_flow"] = std::sqrt((te + te * ti_te) * 1.602176634e-19
		/ (2.014 * 1.66053906660e-27));

	// slab_b_angle is the angle to the wall
	inpts["slab_b_angle"] = 90.0 - theta;

	flan(inpts);

	MPI_Finalize();
	return 0;
}
