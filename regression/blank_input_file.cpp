#include "flan.h"
#include "mpi.h"


Inputs create_inputs()
{
	Inputs inpts;
	
	// Small selection of options to get started with. See documentation for 
	// list of all options: https://shawnzamperini.github.io/flan/input-options/
	
	// Meta options
	inpts["use_gpu"] = "off";  // off, cuda
	inpts["case_name"] = "flan_case_name";
	
	// Impurity options
	inpts["imp_num"] = 100000000;
	inpts["imp_time_step"] = 1e-8;
	inpts["imp_collisions"] = "nanbu";  // off, nanbu

	//   |  Z |  Mass  |
	//He |  2 |   4.00 |
	//inpts["imp_atom_num"] = 2; inpts["imp_mass_amu"] = 4.00;
	//Li |  3 |   6.94 |
	//inpts["imp_atom_num"] = 3; inpts["imp_mass_amu"] = 6.94;
	//B  |  5 |  10.81 |
	//inpts["imp_atom_num"] = 5; inpts["imp_mass_amu"] = 10.81;
	//C  |  6 |  12.01 |
	//inpts["imp_atom_num"] = 6; inpts["imp_mass_amu"] = 12.01;
	//Ne | 10 |  20.18 |
	//inpts["imp_atom_num"] = 10; inpts["imp_mass_amu"] = 20.18;
	//Si | 14 |  20.85 |
	//inpts["imp_atom_num"] = 14; inpts["imp_mass_amu"] = 20.85;
	//Fe | 26 |  55.85 |
	//inpts["imp_atom_num"] = 26; inpts["imp_mass_amu"] = 55.85;
	//Mo | 42 |  95.95 |
	//inpts["imp_atom_num"] = 42; inpts["imp_mass_amu"] = 95.95;
	//Xe | 54 | 131.29 |
	//inpts["imp_atom_num"] = 54; inpts["imp_mass_amu"] = 131.29;
	//W  | 74 | 183.84 |
	inpts["imp_atom_num"] = 74; inpts["imp_mass_amu"] = 183.84;

	// Atomic physics
	inpts["imp_iz_recomb"] = "on";
	inpts["openadas_root"] = "/global/homes/z/zamp/flansoft/openadas";
	inpts["openadas_year"] = 89;  // 50 for tungsten, 89 for everything else


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
