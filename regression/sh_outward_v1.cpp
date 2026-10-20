#include <iostream>
#include <cmath>

#include "flan.h"
#include "flan_types.h"
#include "mpi.h"


Inputs create_inputs()
{
	Inputs inpts;
	
	// Int, double and string input options go here
	inpts["case_name"] = "sh_outward_v1";
	inpts["imp_num"] = 100;
	inpts["gkyl_dir"] = "/home/zamp/gkyldir/d3d-167196-v6-gpu";
	inpts["gkyl_casename"] = "d3d-167196-v6-gpu";
	inpts["gkyl_frame_start"] = 600;
	inpts["gkyl_frame_end"] = 999;
	inpts["gkyl_elec_name"] = "elc";
	inpts["gkyl_ion_name"] = "ion";
	inpts["imp_mass_amu"] = 183.84;
	inpts["imp_init_charge"] = 3;
	inpts["imp_time_step_opt"] = "constant";
	inpts["imp_time_step"] = 1e-9;
	inpts["imp_time_step_min"] = 1e-12;
	inpts["imp_collisions"] = "off";
	inpts["imp_vel_stats"] = "on";
	inpts["imp_xmin"] = 2.315;
	inpts["imp_xmax"] = 2.315;
	inpts["imp_xbound_buffer"] = 0.000;
	inpts["imp_ystart_opt"] = "range";
	inpts["imp_ystart_val"] = -0.10;
	inpts["imp_zstart_val"] = 0.01;
	inpts["imp_iz_recomb"] = "on";
	inpts["openadas_root"] = "/home/zamp/flandir/openadas";
	inpts["openadas_year"] = 50;
	inpts["imp_var_reduct"] = "on";
	inpts["imp_var_reduct_freq"] = 0.05;
	inpts["imp_var_reduct_min_weight"] = 0.9999;

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
