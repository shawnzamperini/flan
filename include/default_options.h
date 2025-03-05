/**
* @file default_options.cpp
* @brief This file contains all the input values and their defaults
*
* The input values and functions are defined in this file along with their
* default values and definitions. They are then overwritten by the
* corresponding values and definitions within the input file that is compiled
* in for a given simulation.
*/
#include <string>

namespace Options
{
	struct Defaults
	{

		// Options related to the Gkeyll simulation
		std::string gkyl_dir {"undefined"};
		std::string gkyl_casename {"undefined"};
		int gkyl_frame_start {0};
		int gkyl_frame_end {1};
		std::string gkyl_elec_name {"elc"};
		std::string gkyl_ion_name {"ion"};
		double gkyl_elec_mass_amu {0.000548};
		double gkyl_ion_mass_amu {2.014};
		//gkyl_file_type,

		// Impurity characteristics
		int imp_atom_num {74};
		double imp_mass_amu {183.84};
		int imp_init_charge {1};

		// Impurity transport simulation settings
		int imp_num {1};
		double imp_xmin {0.0};
		double imp_xmax {0.0};
		std::string imp_ystart_opt {"single_value"};
		double imp_ystart_val {0.0};
		std::string imp_zstart_opt {"single_value"};
		double imp_zstart_val {0.0};
		std::string imp_collisions {"off"};
		std::string imp_var_reduct {"off"};
		double imp_var_reduct_freq {0.1};
		double imp_var_reduct_min_weight {0.1};
		std::string imp_time_step_opt {"variable"};	
		double imp_time_step {1e-7};
		double imp_time_step_min {1e-12};
		double imp_source_scale_fact {1.0};
		std::string imp_vel_stats {"on"};
		double imp_xbound_buffer {0.0};
		std::string imp_iz_recomb {"on"};

		// OpenADAS options
		std::string openadas_root {"undefined"};
		int openadas_year {50};

		// Functions to map between computational and physical (Cartesian) 
		// coordinates. These are defined as if the computational coordinates are
		// already in Cartesian, but the pointer to each can be reseated to one
		// from the simulation input file if necessary.
	}
}
