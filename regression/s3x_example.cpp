#include "flan.h"
#include "mpi.h"


// Example of running Flan on a SOLEDGE3X background: DIII-D #167196, with
// tungsten starting in the near SOL at the outer midplane. See the
// SOLEDGE3X page of the documentation for how SOLEDGE3X backgrounds work.
Inputs create_inputs()
{
	Inputs inpts;

	// Meta options. Multi-zone backgrounds only run on the CPU for now.
	inpts["use_gpu"] = "off";
	inpts["case_name"] = "s3x_example";

	// SOLEDGE3X background. s3x_dir holds the mesh, the reference
	// (normalization) file and the plasma files plasma_00000.h5, ...
	inpts["bkg_source"] = "soledge3x";
	inpts["s3x_dir"] = "/mnt/c/Users/Shawn Zamperini/Documents/"
		"perlmutter_staging/diiid_167196_3d";  // Change to your directory
	inpts["s3x_mesh_file"] = "mesh_S3X.h5";
	inpts["s3x_ref_file"] = "refParam_raptorX.h5";
	inpts["s3x_frame_start"] = 0;
	inpts["s3x_frame_end"] = 1;
	inpts["s3x_chi_absorb"] = "on";

	// The two frames are only 0.5 us apart, so loop over them in time to
	// treat the background as steady
	inpts["tbound_type"] = "periodic";

	// Impurity options
	inpts["imp_num"] = 1000;
	inpts["imp_time_step"] = 1e-9;
	inpts["imp_max_time"] = 1e-3;
	inpts["imp_collisions"] = "nanbu";
	inpts["imp_atom_num"] = 74;
	inpts["imp_mass_amu"] = 183.84;
	inpts["imp_init_charge"] = 1;
	inpts["imp_temp_start_opt"] = "single_value";
	inpts["imp_temp_start_val"] = 5.0;

	// Start in zone 7 (the near SOL), 5.5 cells out in psi, at the theta
	// cell closest to the outer midplane (R = 2.27 m, Z = 0.06 m), and at a
	// random toroidal angle and time. x, y and z are the zone's psi index,
	// theta index and toroidal angle.
	inpts["imp_zone_start_val"] = 7;
	inpts["imp_tstart_opt"] = "full_range";
	inpts["imp_xstart_opt"] = "single_value";
	inpts["imp_xstart_val"] = 5.5;
	inpts["imp_ystart_opt"] = "single_value";
	inpts["imp_ystart_val"] = 473.5;
	inpts["imp_zstart_opt"] = "full_range";

	// Sputtering of the wall by the background plasma. With "sputter", the
	// particles are the W atoms D+ sputters off the wall instead: they start
	// neutral where the sputtered source puts them, the start options above
	// aren't used, and the density is scaled by the source's total rate.
	// The source is also saved in the sputtering group of the output file.
	inpts["imp_source_opt"] = "manual";  // "sputter" to start from the wall
	inpts["sput_table_file"] = "/path/to/flan/data/sputtering/"
		"D_on_W.h5";  // Change to your directory
	inpts["sput_flux_opt"] = "parallel_exb";  // Parallel flow + E x B drift
	inpts["sput_bohm_floor"] = "off";  // "on": parallel flow >= sound speed
	inpts["sput_angle_val"] = 60.0;  // Ion angle of incidence (deg)
	inpts["sput_max_wall_dist"] = 2.0;  // Cell widths to the wall outline

	// Atomic physics
	inpts["imp_iz_recomb"] = "on";
	inpts["openadas_root"] = "/path/to/openadas";  // Change to your directory
	inpts["openadas_year"] = 50;  // 50 for tungsten, 89 for everything else

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
