#include "flan.h"
#include "mpi.h"


// Example of running Flan on a HERMES (JAX-Hermes3) background: a
// lower single-null turbulence run, with tungsten starting just outside
// the separatrix at the outer midplane. See the HERMES page of the
// documentation for how HERMES backgrounds work.
Inputs create_inputs()
{
	Inputs inpts;

	// Meta options. Multi-zone backgrounds only run on the CPU for now.
	inpts["use_gpu"] = "off";
	inpts["case_name"] = "hermes_example";

	// HERMES background: the run's hypnotoad grid and its directory of
	// snapshots (snap*.npz). Load the last 50 snapshots.
	inpts["bkg_source"] = "hermes";
	inpts["hermes_format"] = "jax";
	inpts["hermes_grid_file"] = "/path/to/run/mesh/grid.nc";  // Change
	inpts["hermes_snap_dir"] = "/path/to/run/snapshots";  // Change
	inpts["hermes_snap_start"] = -50;
	inpts["hermes_snap_end"] = -1;
	inpts["hermes_snap_stride"] = 1;

	// These must be the run's, from the case's input JSON: pars.norm
	// (Nnorm, Tnorm, Bnorm), pars.average_atomic_mass, pars.zperiod and
	// pars.kappa_coefficient_e and _i
	inpts["hermes_nnorm"] = 1e19;
	inpts["hermes_tnorm"] = 50.0;
	inpts["hermes_bnorm"] = 1.0;
	inpts["hermes_ion_mass_amu"] = 2.0;
	inpts["hermes_zperiod"] = 5;
	inpts["hermes_kappa_coef_e"] = 2.23445742854949;  // 3.16 / sqrt(2)
	inpts["hermes_kappa_coef_i"] = 3.9;

	// The snapshots hold no heat fluxes, so the thermal force uses ones
	// calculated from the conduction. The plasma is followed along the
	// field line ("off" is only for comparison).
	inpts["hermes_heat_flux"] = "conduction";
	inpts["hermes_field_aligned"] = "on";

	// Snapshots are a short stretch of a turbulent run, so loop over them
	// in time and treat the turbulence as statistically steady
	inpts["tbound_type"] = "periodic";

	// Impurity options
	inpts["imp_num"] = 1000;
	inpts["imp_time_step"] = 1e-9;
	inpts["imp_max_time"] = 1e-3;
	inpts["imp_collisions"] = "nanbu";
	inpts["imp_thermal_force"] = "on";
	inpts["imp_atom_num"] = 74;
	inpts["imp_mass_amu"] = 183.84;
	inpts["imp_init_charge"] = 1;
	inpts["imp_temp_start_opt"] = "single_value";
	inpts["imp_temp_start_val"] = 5.0;

	// Start in zone 5 (the main SOL) at the outer midplane, at a random
	// toroidal angle and time. x and y are the grid's global cell indices
	// (x counting the two guard cells, so the separatrix is at x =
	// ixseps1) and z is the toroidal angle in the wedge. For the grid this
	// was set up for, ixseps1 = 18 and Z = 0 at the outer midplane is at
	// y = 18.25: check these on your grid (Rxy, Zxy).
	inpts["imp_zone_start_val"] = 5;
	inpts["imp_tstart_opt"] = "full_range";
	inpts["imp_xstart_opt"] = "single_value";
	inpts["imp_xstart_val"] = 20.5;
	inpts["imp_ystart_opt"] = "single_value";
	inpts["imp_ystart_val"] = 18.25;
	inpts["imp_zstart_opt"] = "full_range";
	inpts["imp_source_opt"] = "manual";

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
