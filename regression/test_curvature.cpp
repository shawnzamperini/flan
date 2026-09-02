#include <iostream>
#include <cmath>
#include <tuple>

#include "flan.h"
#include "flan_types.h"


// mapc2p (slated for depreciation)
std::tuple<double, double, double> mapc2p(double xc, double yc, double zc)
{
	std::cout << "Don't use me!\n";
	return {0, 0 ,0};
}

Inputs create_inputs()
{
	Inputs inpts;
	
	// Setup the ExB drift test
	inpts["bkg_source"] = "test";
	inpts["case_name"] = "curvature";
	inpts["test_opt"] = "curvature";
	inpts["save_track"] = "on";  // Needed for validating tests

	// Impurity and transport options
	inpts["imp_num"] = 1;  // Only 1 particle for drift tests

	/*
	   |  Z |  Mass  |
	He |  2 |   4.00 |
	Li |  3 |   6.94 |
	B  |  5 |  10.81 |
	C  |  6 |  12.01 |
	Ne | 10 |  20.18 |
	Fe | 26 |  55.85 |
	Mo | 42 |  95.95 |
	W  | 74 | 183.84 |
	*/
	inpts["imp_atom_num"] = 74;
	inpts["imp_mass_amu"] = 183.84;
	inpts["imp_init_charge"] = 5;
	inpts["imp_time_step_opt"] = "constant";
	inpts["imp_time_step"] = 1e-8;
	inpts["imp_time_step_min"] = 1e-12;
	inpts["imp_collisions"] = "off";

	// Starting location options. We want to start particles at the same time
	// and place so we can track their motion as one single particle. This is
	// good because if there was numerical noise we'd see it in the videos if
	// the particles didn't stay mostly together.
	inpts["imp_tstart_opt"] = "single_value";
	inpts["imp_tstart_val"] = 0.0;
	inpts["imp_xstart_opt"] = "single_value";
	inpts["imp_xstart_val"] = 2.01;  // x range from 2.00 - 2.05
	inpts["imp_ystart_opt"] = "single_value";
	inpts["imp_ystart_val"] = 0.0;
	inpts["imp_zstart_opt"] = "single_value";
	inpts["imp_zstart_val"] = 0.0;  // z range from 0 - 2pi

	// These get overwritten in the test, just leaving here to make this note
	inpts["imp_temp_start_opt"] = "single_value";
	inpts["imp_temp_start_val"] = 50.0;

	// Atomic physics, want this off
	inpts["imp_iz_recomb"] = "off";

	// Miscellanous
	inpts["print_interval"] = 10;

	// Pointers to functions go here
	Mapc2p_ptr mapc2p_ptr {&mapc2p};

	// Function pointer input options go here
	inpts["mapc2p"] = mapc2p_ptr;

	return inpts;
}


int main()
{
	// Run Flan
	Inputs inpts {create_inputs()};
	flan(inpts);
}
