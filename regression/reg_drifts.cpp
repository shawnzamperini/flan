/**
* @file reg_drifts.cpp
*
* @brief Regression cases: gyration and the guiding-center drifts
*
* Singly charged tungsten ions (W5+ for the curvature drift) are followed
* without collisions in simple fields, and their motion is compared with
* the textbook results for a single particle. Every particle of a run
* starts at the same place with the same velocity, so they all do the
* same thing and the output's mean velocity is that of one particle.
*
* Usage: reg_drifts <variant>, where <variant> is one of
*
*   gyration      B only. The particle gyrates at the cyclotron frequency
*                 q B / m and keeps its speed.
*   exb           B and a steady E across it. Drift of E / B, across both.
*   polarization  B and an E that grows steadily. Polarization drift of
*                 (m / (q B^2)) dE/dt along E, on top of the ExB drift.
*   gradb         B that gets stronger across the field. Grad-B drift of
*                 m v_perp^2 / (2 q B) times the relative gradient.
*   curvature     Toroidal B (falling as 1 / R) and a particle moving
*                 along it. Curvature drift of m v_par^2 / (q B R),
*                 vertically.
*   gradb_gpu     As gradb, on the GPU.
*
* In the straight slab B is along +X, so for a positive ion: gyration
* turns +Y toward -Z, an E along +Y gives an ExB drift along -Z, and a
* gradient along +Y gives a grad-B drift along +Z. In the bent slab B is
* along the toroidal angle and the curvature drift is along +Z (upward).
* check_drifts.py reads the output, works out the expected motion from
* the inputs saved in it, and compares.
*
* Everything about a case is set in this file: the fields come from the
* slab background (bkg_source = "slab") and the starting velocity from
* imp_v*_start_val. The plasma density and temperatures are there only
* because a background needs them; with collisions and ionization off
* they do nothing.
*
* The output is reg_drifts_<variant>.nc in the working directory.
*/
#include <cmath>
#include <iostream>
#include <string>

#include "flan.h"
#include "flan_types.h"
#include "mpi.h"

/**
* Inputs shared by every variant: a straight slab with B = 1 T along +X
* and no other field, periodic boundaries, and 100 W+ ions starting at
* rest in the middle of the box at t = 0, followed for 100 us with a time
* step of 1e-8 s (1/1200 of a gyration) and no collisions. The variants
* change what they need.
*
* @param variant Name of the variant, used in the case name.
* @return The inputs.
*/
Inputs base_inputs(const std::string& variant)
{
	Inputs inpts;

	inpts["case_name"] = "reg_drifts_" + variant;
	inpts["use_gpu"] = "off";

	// Fields. The plasma values are placeholders (see the top of the
	// file).
	inpts["bkg_source"] = "slab";
	inpts["gkyl_ion_mass_amu"] = 2.014;
	inpts["gkyl_elec_mass_amu"] = 0.000548;
	inpts["slab_ne"] = 1e19;
	inpts["slab_te"] = 10.0;
	inpts["slab_ti"] = 10.0;
	inpts["slab_flow"] = 0.0;
	inpts["slab_heat_flux_i"] = 0.0;
	inpts["slab_heat_flux_e"] = 0.0;
	inpts["slab_bmag"] = 1.0;
	inpts["slab_efield_y"] = 0.0;
	inpts["slab_efield_y_rate"] = 0.0;
	inpts["slab_bgrad_y"] = 0.0;
	inpts["slab_major_radius"] = 0.0;
	inpts["slab_tmax"] = 100e-6;

	// Particles leaving one side of the box come back in the other
	inpts["min_xbound_type"] = "periodic";
	inpts["max_xbound_type"] = "periodic";
	inpts["min_ybound_type"] = "periodic";
	inpts["max_ybound_type"] = "periodic";
	inpts["min_zbound_type"] = "periodic";
	inpts["max_zbound_type"] = "periodic";

	// Tungsten: W+ gyrates in 12 us at 1 T
	inpts["imp_num"] = 100;
	inpts["imp_atom_num"] = 74;
	inpts["imp_mass_amu"] = 183.84;
	inpts["imp_init_charge"] = 1;
	inpts["imp_iz_recomb"] = "off";
	inpts["imp_tstart_opt"] = "single_value";
	inpts["imp_tstart_val"] = 0.0;
	inpts["imp_xstart_opt"] = "single_value";
	inpts["imp_xstart_val"] = 0.025;
	inpts["imp_ystart_opt"] = "single_value";
	inpts["imp_ystart_val"] = 0.0;
	inpts["imp_zstart_opt"] = "single_value";
	inpts["imp_zstart_val"] = 0.0;

	// Starting velocity: no thermal part, so exactly imp_v*_start_val
	inpts["imp_temp_start_opt"] = "single_value";
	inpts["imp_temp_start_val"] = 0.0;
	inpts["imp_vX_start_val"] = 0.0;
	inpts["imp_vY_start_val"] = 0.0;
	inpts["imp_vZ_start_val"] = 0.0;

	inpts["imp_time_step_opt"] = "constant";
	inpts["imp_time_step"] = 1e-8;
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

	if (variant == "gyration")
	{
		// 10 km/s across the field: a gyroradius of 1.9 cm. About four
		// gyrations, so the output's time frames resolve them.
		inpts["imp_vY_start_val"] = 10000.0;
		inpts["slab_tmax"] = 50e-6;
	}
	else if (variant == "exb")
	{
		// 500 V/m across 1 T: a drift of 500 m/s. Starting at rest, the
		// particle also gyrates at 500 m/s about the drifting center.
		inpts["slab_efield_y"] = 500.0;
	}
	else if (variant == "polarization")
	{
		// E grows from 0 by 1e8 V/m per second (100 V/m per us): a
		// polarization drift of 190 m/s along E, and an ExB drift that
		// grows to 10 km/s by the end
		inpts["slab_efield_y_rate"] = 1e8;
	}
	else if (variant == "gradb" || variant == "gradb_gpu")
	{
		// B grows by 2 T/m (a relative gradient of 2 /m) across the
		// field, and the particle gyrates at 10 km/s: a drift of 191
		// m/s. The gyroradius is 4% of the gradient length, so the
		// formula holds to a fraction of a percent. The field isn't
		// periodic in y; the particle stays within 2 cm of y = 0.
		inpts["slab_bgrad_y"] = 2.0;
		inpts["imp_vY_start_val"] = 10000.0;
		inpts["min_ybound_type"] = "absorbing";
		inpts["max_ybound_type"] = "absorbing";
		if (variant == "gradb_gpu") inpts["use_gpu"] = "cuda";
	}
	else if (variant == "curvature")
	{
		// Slab bent to a major radius of 2 m: toroidal field of 1 T at
		// R = 2 m. W5+ starts at R = 2 m, moving along the field at 30
		// km/s: a drift of 172 m/s upward. In 80 us it goes 1.2 rad
		// around, staying inside the quarter torus (0 to pi/2), and 1.4
		// cm up. Coordinates are x = R, y = Z, z = toroidal angle; at
		// the starting angle of 0.05 the toroidal direction is (-sin,
		// cos, 0) of it in X, Y, Z.
		constexpr double v_par {30000.0};
		constexpr double phi_start {0.05};
		inpts["slab_major_radius"] = 2.0;
		inpts["slab_tmax"] = 80e-6;
		inpts["imp_init_charge"] = 5;
		inpts["imp_xstart_val"] = 2.0;
		inpts["imp_ystart_val"] = 0.0;
		inpts["imp_zstart_val"] = phi_start;
		inpts["imp_vX_start_val"] = -v_par * std::sin(phi_start);
		inpts["imp_vY_start_val"] = v_par * std::cos(phi_start);
		inpts["min_xbound_type"] = "absorbing";
		inpts["max_xbound_type"] = "absorbing";
		inpts["min_ybound_type"] = "absorbing";
		inpts["max_ybound_type"] = "absorbing";
		inpts["min_zbound_type"] = "absorbing";
		inpts["max_zbound_type"] = "absorbing";
	}
	else
	{
		std::cerr << "Usage: reg_drifts <variant>, with <variant> one of "
			"gyration, exb, polarization, gradb, curvature, gradb_gpu\n";
		MPI_Finalize();
		return 1;
	}

	flan(inpts);

	MPI_Finalize();
	return 0;
}
