#include <iostream>
#include <cmath>

#include "flan.h"
#include "flan_types.h"
#include "mpi.h"


Inputs create_inputs()
{
	Inputs inpts;
	
	// Setup the ExB drift test
	inpts["bkg_source"] = "test";
	inpts["case_name"] = "gradb";
	inpts["test_opt"] = "gradb";

	// Impurity and transport options
	inpts["imp_num"] = 100;

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
	inpts["imp_xstart_val"] = 0.01;
	inpts["imp_ystart_opt"] = "single_value";
	inpts["imp_ystart_val"] = 0.0;
	inpts["imp_zstart_opt"] = "single_value";
	inpts["imp_zstart_val"] = 0.0;

	// These get overwritten in the test, just leaving here to make this note
	inpts["imp_temp_start_opt"] = "single_value";
	inpts["imp_temp_start_val"] = 50.0;

	// Atomic physics, want this off
	inpts["imp_iz_recomb"] = "off";

	// Miscellanous
	inpts["print_interval"] = 10;

	return inpts;
}


int main(int argc, char* argv[])
{
	// Initialize MPI
	MPI_Init(&argc, &argv);

	// Bundle up inputs and run Flan
	Inputs inpts {create_inputs()};
	flan(inpts);

	// Finalize MPI
	MPI_Finalize();
}
