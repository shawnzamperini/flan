/**
* @file options.cpp
* @brief Function definitoins for Options class.
*/
#include <iostream>
#include <string>
#include <tuple>
#include <initializer_list>

#include "options.h"
#include "options_device.h"


namespace Options
{

	/**
	* Constructs an Options object with every option at its default value.
	*
	* The defaults are set where each member variable is declared in
	* options.h. The control integers are then set to match the default
	* string options by calling initialize_control_ints.
	*/
	Options::Options()
	{
		// Set internal control integers based on (at this point, default)
		// options. They'll get changed later if the option changes.
		initialize_control_ints();
	}
	
	/**
	* Sets every control integer to match its string option.
	*
	* Control integers are integer versions of string options that are
	* faster
	* to check in the code than comparing strings. The logic that assigns
	* each one lives in the option's setter, so this calls every such setter
	* with the option's current value. Any new option with a control integer
	* must have its setter called here too. Otherwise its control integer
	* stays at -1, which get_control_int reports as an error.
	*/
	void Options::initialize_control_ints()
	{
		// We are setting each value to the value it is already at, why?
		// Because the logic for assigning the respective control integers
		// lives within each set_ function. This is not the best, since if
		// one were to add an additional control integer for a new option and
		// then forgot to call the corresponding set_ function within this
		// function, then they may introduce a silent bug! This can be silent
		// because if the user does not have that option in their input file,
		// the default string value will be used but there is no guarantee that
		// the default control integer will be the correct value! Luckily,
		// this bug is avoided by accessing the control integers via
		// get_control_int, which will tell you if a control integer hasn't
		// been correctly initalized (by not including it here).
		set_use_gpu(m_use_gpu);
		set_bkg_source(m_bkg_source);
		set_test_opt(m_test_opt);
		set_save_track(m_save_track);
		set_imp_tstart_opt(m_imp_tstart_opt);
		set_imp_xstart_opt(m_imp_xstart_opt);
		set_imp_ystart_opt(m_imp_ystart_opt);
		set_imp_zstart_opt(m_imp_zstart_opt);
		set_imp_temp_start_opt(m_imp_temp_start_opt);
		set_imp_collisions(m_imp_collisions);
		set_var_red_split(m_var_red_split);
		set_var_red_import(m_var_red_import);
		set_var_red_rusrol(m_var_red_rusrol);
		set_imp_time_step_opt(m_imp_time_step_opt);
		set_imp_iz_recomb(m_imp_iz_recomb);
		set_tbound_type(m_tbound_type);
		set_min_xbound_type(m_min_xbound_type);
		set_max_xbound_type(m_max_xbound_type);
		set_min_ybound_type(m_min_ybound_type);
		set_max_ybound_type(m_max_ybound_type);
		set_min_zbound_type(m_min_zbound_type);
		set_max_zbound_type(m_max_zbound_type);
		set_calc_grad_elec(m_calc_grad_elec);
	}

	// Setter definitions. They're all the same, so just throwing together
	// as a mass of code. Can organize if necessary, but it would probably
	// just take up three times as much space.

	/**
	* Sets the case_name option.
	*
	* Name of the simulation case. The results are saved to <case_name>.nc.
	*
	* @param case_name New value of the case_name option.
	*/
	void Options::set_case_name(std::string case_name) 
		{m_case_name = case_name;}

	/**
	* Sets the use_gpu option.
	*
	* Whether to run the simulation on GPUs.
	*
	* Valid values are "off" and "cuda". An invalid value prints an error
	* and leaves the option unchanged. Also sets the control integer
	* returned by use_gpu_int() to 0 for "off", 1 for "cuda". An invalid
	* value leaves the control integer unchanged.
	*
	* @param use_gpu New value of the use_gpu option.
	*/
	void Options::set_use_gpu(std::string use_gpu) 
	{
		if (check_input<std::string>("use_gpu", use_gpu,
			{"off", "cuda"}))
		{
			m_use_gpu = use_gpu;
		}

		// Assign control integers
		if (use_gpu == "off") m_use_gpu_int = 0;
		else if (use_gpu == "cuda") m_use_gpu_int = 1;
	}

	/**
	* Sets the seed option.
	*
	* Base seed for the random number generators. A negative value seeds
	* them from the clock so each run is different. The MPI rank is added to
	* the seed so each rank gets a different sequence.
	*
	* @param seed New value of the seed option.
	*/
	void Options::set_seed(int seed) 
		{m_seed = seed;}

	/**
	* Sets the slot_cap option.
	*
	* Maximum number of particle slots, i.e., the most particles followed at
	* once. The number of slots used is the smaller of this and imp_num.
	*
	* @param slot_cap New value of the slot_cap option.
	*/
	void Options::set_slot_cap(int slot_cap) 
		{m_slot_cap = slot_cap;}

	/**
	* Sets the bkg_source option.
	*
	* Source of the background plasma, either a built-in test background or
	* a Gkeyll simulation.
	*
	* Valid values are "test" and "gkeyll". An invalid value prints an error
	* and leaves the option unchanged. Also sets the control integer
	* returned by bkg_source_int() to 0 for "test", 1 for "gkeyll". An
	* invalid value leaves the control integer unchanged.
	*
	* @param bkg_source New value of the bkg_source option.
	*/
	void Options::set_bkg_source(std::string bkg_source) 
	{
		if (check_input<std::string>("bkg_source", bkg_source,
			{"test", "gkeyll"}))
		{
			m_bkg_source = bkg_source;
		}

		// Assign control integers
		if (bkg_source == "test") m_bkg_source_int = 0;
		else if (bkg_source == "gkeyll") m_bkg_source_int = 1;
	}

	/**
	* Sets the test_opt option.
	*
	* Which built-in test case to run when using the test background.
	*
	* Valid values are "gyrate", "exb", "gradb", "polarization", "curvature"
	* and "friction_force". An invalid value prints an error and leaves the
	* option unchanged. Also sets the control integer returned by
	* test_opt_int() to 0 for "gyrate", 1 for "exb", 2 for "gradb", 3 for
	* "polarization", 4 for "curvature", 5 for "friction_force". An invalid
	* value leaves the control integer unchanged.
	*
	* @param test_opt New value of the test_opt option.
	*/
	void Options::set_test_opt(std::string test_opt) 
	{
		if (check_input<std::string>("test_opt", test_opt,
			{"gyrate", "exb", "gradb", "polarization", "curvature", 
			"friction_force"}))
		{
			m_test_opt = test_opt;
		}

		// Assign control integers
		if (test_opt == "gyrate") m_test_opt_int = 0;
		else if (test_opt == "exb") m_test_opt_int = 1;
		else if (test_opt == "gradb") m_test_opt_int = 2;
		else if (test_opt == "polarization") m_test_opt_int = 3;
		else if (test_opt == "curvature") m_test_opt_int = 4;
		else if (test_opt == "friction_force") m_test_opt_int = 5;
	}

	/**
	* Sets the save_track option.
	*
	* Whether to save particle tracks to the output file.
	*
	* Valid values are "off" and "on". An invalid value prints an error and
	* leaves the option unchanged. Also sets the control integer returned by
	* save_track_int() to 0 for "off", 1 for "on". An invalid value leaves
	* the control integer unchanged.
	*
	* @param save_track New value of the save_track option.
	*/
	void Options::set_save_track(std::string save_track) 
	{
		if (check_input<std::string>("save_track", save_track,
			{"off", "on"}))
		{
			m_save_track = save_track;
		}

		// Assign control integers
		if (save_track == "off") m_save_track_int = 0;
		else if (save_track == "on") m_save_track_int = 1;
	}
	
	/**
	* Sets the gkyl_dir option.
	*
	* Directory containing the Gkeyll simulation files.
	*
	* @param gkyl_dir New value of the gkyl_dir option.
	*/
	void Options::set_gkyl_dir(std::string gkyl_dir) 
		{m_gkyl_dir = gkyl_dir;}
	
	/**
	* Sets the gkyl_casename option.
	*
	* Name of the Gkeyll case, i.e., the prefix of its output files.
	*
	* @param gkyl_casename New value of the gkyl_casename option.
	*/
	void Options::set_gkyl_casename(std::string gkyl_casename) 
		{m_gkyl_casename = gkyl_casename;}
	
	/**
	* Sets the gkyl_frame_start option.
	*
	* First Gkeyll frame to load.
	*
	* @param gkyl_frame_start New value of the gkyl_frame_start option.
	*/
	void Options::set_gkyl_frame_start(int gkyl_frame_start) 
		{m_gkyl_frame_start = gkyl_frame_start;}
	
	/**
	* Sets the gkyl_frame_end option.
	*
	* Last Gkeyll frame to load. The range from gkyl_frame_start to
	* gkyl_frame_end is inclusive.
	*
	* @param gkyl_frame_end New value of the gkyl_frame_end option.
	*/
	void Options::set_gkyl_frame_end(int gkyl_frame_end) 
		{m_gkyl_frame_end = gkyl_frame_end;}

	/**
	* Sets the gkyl_elec_name option.
	*
	* Name of the electron species in the Gkeyll simulation.
	*
	* @param gkyl_elec_name New value of the gkyl_elec_name option.
	*/
	void Options::set_gkyl_elec_name(std::string gkyl_elec_name) 
		{m_gkyl_elec_name = gkyl_elec_name;}
	
	/**
	* Sets the gkyl_ion_name option.
	*
	* Name of the main ion species in the Gkeyll simulation.
	*
	* @param gkyl_ion_name New value of the gkyl_ion_name option.
	*/
	void Options::set_gkyl_ion_name(std::string gkyl_ion_name) 
		{m_gkyl_ion_name = gkyl_ion_name;}
	
	/**
	* Sets the gkyl_elec_mass_amu option.
	*
	* Electron mass used in the Gkeyll simulation (amu).
	*
	* @param gkyl_elec_mass_amu New value of the gkyl_elec_mass_amu option.
	*/
	void Options::set_gkyl_elec_mass_amu(double gkyl_elec_mass_amu) 
		{m_gkyl_elec_mass_amu = gkyl_elec_mass_amu;}
	
	/**
	* Sets the gkyl_ion_mass_amu option.
	*
	* Main ion mass used in the Gkeyll simulation (amu).
	*
	* @param gkyl_ion_mass_amu New value of the gkyl_ion_mass_amu option.
	*/
	void Options::set_gkyl_ion_mass_amu(double gkyl_ion_mass_amu) 
		{m_gkyl_ion_mass_amu = gkyl_ion_mass_amu;}

	/**
	* Sets the gkyl_file_type option.
	*
	* File format of the Gkeyll output.
	*
	* Valid values are "binary". An invalid value prints an error and leaves
	* the option unchanged.
	*
	* @param gkyl_file_type New value of the gkyl_file_type option.
	*/
	void Options::set_gkyl_file_type(std::string gkyl_file_type) 
	{
		if (check_input<std::string>("gkyl_file_type", gkyl_file_type,
			{"binary"}))
		{
			m_gkyl_file_type = gkyl_file_type;
		}
	}

	/**
	* Sets the gkyl_moment_type option.
	*
	* Type of moment files in the Gkeyll output. Collisions need
	* "bimaxwellian", which provides the ion parallel flow and perpendicular
	* velocity used to calculate the background flow.
	*
	* Valid values are "bimaxwellian", "maxwellian" and "m0m1m2". An invalid
	* value prints an error and leaves the option unchanged.
	*
	* @param gkyl_moment_type New value of the gkyl_moment_type option.
	*/
	void Options::set_gkyl_moment_type(std::string gkyl_moment_type) 
	{
		if (check_input<std::string>("gkyl_moment_type", gkyl_moment_type,
			{"bimaxwellian", "maxwellian", "m0m1m2"}))
		{
			m_gkyl_moment_type = gkyl_moment_type;
		}
	}

	/**
	* Sets the calc_grad_elec option.
	*
	* Whether to calculate the electric field gradient. This is not
	* implemented yet, so turning it on only prints a warning.
	*
	* Valid values are "off" and "on". An invalid value prints an error and
	* leaves the option unchanged. Also sets the control integer returned by
	* calc_grad_elec_int() to 0 for "off", 1 for "on". An invalid value
	* leaves the control integer unchanged.
	*
	* @param calc_grad_elec New value of the calc_grad_elec option.
	*/
	void Options::set_calc_grad_elec(std::string calc_grad_elec) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("calc_grad_elec", calc_grad_elec, 
				{"off", "on"}))
			{
				m_calc_grad_elec = calc_grad_elec;
			}

			// Assign control integers
			if (calc_grad_elec == "off") m_calc_grad_elec_int = 0;
			else if (calc_grad_elec == "on") m_calc_grad_elec_int = 1;
		}
	
	/**
	* Sets the min_ne option.
	*
	* Lower limit on the electron density (m^-3).
	*
	* @param min_ne New value of the min_ne option.
	*/
	void Options::set_min_ne(double min_ne) 
		{m_min_ne = min_ne;}

	/**
	* Sets the min_te option.
	*
	* Lower limit on the electron temperature (eV).
	*
	* @param min_te New value of the min_te option.
	*/
	void Options::set_min_te(double min_te) 
		{m_min_te = min_te;}

	/**
	* Sets the min_ti option.
	*
	* Lower limit on the ion temperature (eV), applied when calculating
	* collisions with the background ions.
	*
	* @param min_ti New value of the min_ti option.
	*/
	void Options::set_min_ti(double min_ti) 
		{m_min_ti = min_ti;}

	/**
	* Sets the tbound_type option.
	*
	* Boundary condition at the start and end of the background time range.
	*
	* Valid values are "absorbing" and "periodic". An invalid value prints
	* an error and leaves the option unchanged. Also sets the control
	* integer returned by tbound_type_int() to 0 for "absorbing", 1 for
	* "periodic". An invalid value leaves the control integer unchanged.
	*
	* @param tbound_type New value of the tbound_type option.
	*/
	void Options::set_tbound_type(std::string tbound_type) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("tbound_type", tbound_type,
			{"absorbing", "periodic"}))
		{
			m_tbound_type = tbound_type;
		}

		// Assign control integers
		if (tbound_type == "absorbing") m_tbound_type_int = 0;
		else if (tbound_type == "periodic") m_tbound_type_int = 1;
	}

	/**
	* Sets the min_xbound_type option.
	*
	* Boundary condition at the minimum x boundary.
	*
	* Valid values are "absorbing", "periodic" and "core". An invalid value
	* prints an error and leaves the option unchanged. Also sets the control
	* integer returned by min_xbound_type_int() to 0 for "absorbing", 1 for
	* "periodic", 2 for "core". An invalid value leaves the control integer
	* unchanged.
	*
	* @param min_xbound_type New value of the min_xbound_type option.
	*/
	void Options::set_min_xbound_type(std::string min_xbound_type) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("min_xbound_type", min_xbound_type,
			{"absorbing", "periodic", "core"}))
		{
			m_min_xbound_type = min_xbound_type;
		}

		// Assign control integers
		if (min_xbound_type == "absorbing") m_min_xbound_type_int = 0;
		else if (min_xbound_type == "periodic") m_min_xbound_type_int = 1;
		else if (min_xbound_type == "core") m_min_xbound_type_int = 2;
	}

	/**
	* Sets the max_xbound_type option.
	*
	* Boundary condition at the maximum x boundary.
	*
	* Valid values are "absorbing", "periodic" and "core". An invalid value
	* prints an error and leaves the option unchanged. Also sets the control
	* integer returned by max_xbound_type_int() to 0 for "absorbing", 1 for
	* "periodic", 2 for "core". An invalid value leaves the control integer
	* unchanged.
	*
	* @param max_xbound_type New value of the max_xbound_type option.
	*/
	void Options::set_max_xbound_type(std::string max_xbound_type) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("max_xbound_type", max_xbound_type,
			{"absorbing", "periodic", "core"}))
		{
			m_max_xbound_type = max_xbound_type;
		}

		// Assign control integers
		if (max_xbound_type == "absorbing") m_max_xbound_type_int = 0;
		else if (max_xbound_type == "periodic") m_max_xbound_type_int = 1;
		else if (max_xbound_type == "core") m_max_xbound_type_int = 2;
	}

	/**
	* Sets the min_ybound_type option.
	*
	* Boundary condition at the minimum y boundary.
	*
	* Valid values are "absorbing", "periodic" and "core". An invalid value
	* prints an error and leaves the option unchanged. Also sets the control
	* integer returned by min_ybound_type_int() to 0 for "absorbing", 1 for
	* "periodic", 2 for "core". An invalid value leaves the control integer
	* unchanged.
	*
	* @param min_ybound_type New value of the min_ybound_type option.
	*/
	void Options::set_min_ybound_type(std::string min_ybound_type) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("min_ybound_type", min_ybound_type,
			{"absorbing", "periodic", "core"}))
		{
			m_min_ybound_type = min_ybound_type;
		}

		// Assign control integers
		if (min_ybound_type == "absorbing") m_min_ybound_type_int = 0;
		else if (min_ybound_type == "periodic") m_min_ybound_type_int = 1;
		else if (min_ybound_type == "core") m_min_ybound_type_int = 2;
	}

	/**
	* Sets the max_ybound_type option.
	*
	* Boundary condition at the maximum y boundary.
	*
	* Valid values are "absorbing", "periodic" and "core". An invalid value
	* prints an error and leaves the option unchanged. Also sets the control
	* integer returned by max_ybound_type_int() to 0 for "absorbing", 1 for
	* "periodic", 2 for "core". An invalid value leaves the control integer
	* unchanged.
	*
	* @param max_ybound_type New value of the max_ybound_type option.
	*/
	void Options::set_max_ybound_type(std::string max_ybound_type) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("max_ybound_type", max_ybound_type,
			{"absorbing", "periodic", "core"}))
		{
			m_max_ybound_type = max_ybound_type;
		}

		// Assign control integers
		if (max_ybound_type == "absorbing") m_max_ybound_type_int = 0;
		else if (max_ybound_type == "periodic") m_max_ybound_type_int = 1;
		else if (max_ybound_type == "core") m_max_ybound_type_int = 2;
	}

	/**
	* Sets the min_zbound_type option.
	*
	* Boundary condition at the minimum z boundary.
	*
	* Valid values are "absorbing", "periodic" and "core". An invalid value
	* prints an error and leaves the option unchanged. Also sets the control
	* integer returned by min_zbound_type_int() to 0 for "absorbing", 1 for
	* "periodic", 2 for "core". An invalid value leaves the control integer
	* unchanged.
	*
	* @param min_zbound_type New value of the min_zbound_type option.
	*/
	void Options::set_min_zbound_type(std::string min_zbound_type) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("min_zbound_type", min_zbound_type,
			{"absorbing", "periodic", "core"}))
		{
			m_min_zbound_type = min_zbound_type;
		}

		// Assign control integers
		if (min_zbound_type == "absorbing") m_min_zbound_type_int = 0;
		else if (min_zbound_type == "periodic") m_min_zbound_type_int = 1;
		else if (min_zbound_type == "core") m_min_zbound_type_int = 2;
	}

	/**
	* Sets the max_zbound_type option.
	*
	* Boundary condition at the maximum z boundary.
	*
	* Valid values are "absorbing", "periodic" and "core". An invalid value
	* prints an error and leaves the option unchanged. Also sets the control
	* integer returned by max_zbound_type_int() to 0 for "absorbing", 1 for
	* "periodic", 2 for "core". An invalid value leaves the control integer
	* unchanged.
	*
	* @param max_zbound_type New value of the max_zbound_type option.
	*/
	void Options::set_max_zbound_type(std::string max_zbound_type) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("max_zbound_type", max_zbound_type,
			{"absorbing", "periodic", "core"}))
		{
			m_max_zbound_type = max_zbound_type;
		}

		// Assign control integers
		if (max_zbound_type == "absorbing") m_max_zbound_type_int = 0;
		else if (max_zbound_type == "periodic") m_max_zbound_type_int = 1;
		else if (max_zbound_type == "core") m_max_zbound_type_int = 2;
	}

	/**
	* Sets the lcfs_x option.
	*
	* x coordinate of the last closed flux surface.
	*
	* @param lcfs_x New value of the lcfs_x option.
	*/
	void Options::set_lcfs_x(double lcfs_x) 
		{m_lcfs_x = lcfs_x;}

	/**
	* Sets the lcfs_y option.
	*
	* y coordinate of the last closed flux surface.
	*
	* @param lcfs_y New value of the lcfs_y option.
	*/
	void Options::set_lcfs_y(double lcfs_y) 
		{m_lcfs_y = lcfs_y;}

	/**
	* Sets the lcfs_z option.
	*
	* z coordinate of the last closed flux surface.
	*
	* @param lcfs_z New value of the lcfs_z option.
	*/
	void Options::set_lcfs_z(double lcfs_z) 
		{m_lcfs_z = lcfs_z;}

	/**
	* Sets the sep_x_bc_xp_z1 option.
	*
	* First z value for the separatrix x boundary condition near the
	* X-point.
	*
	* @param sep_x_bc_xp_z1 New value of the sep_x_bc_xp_z1 option.
	*/
	void Options::set_sep_x_bc_xp_z1(double sep_x_bc_xp_z1) 
		{m_sep_x_bc_xp_z1 = sep_x_bc_xp_z1;}

	/**
	* Sets the sep_x_bc_xp_z2 option.
	*
	* Second z value for the separatrix x boundary condition near the
	* X-point.
	*
	* @param sep_x_bc_xp_z2 New value of the sep_x_bc_xp_z2 option.
	*/
	void Options::set_sep_x_bc_xp_z2(double sep_x_bc_xp_z2) 
		{m_sep_x_bc_xp_z2 = sep_x_bc_xp_z2;}

	/**
	* Sets the imp_atom_num option.
	*
	* Atomic number of the impurity.
	*
	* @param imp_atom_num New value of the imp_atom_num option.
	*/
	void Options::set_imp_atom_num(int imp_atom_num) 
		{m_imp_atom_num = imp_atom_num;}

	/**
	* Sets the imp_mass_amu option.
	*
	* Mass of the impurity (amu).
	*
	* @param imp_mass_amu New value of the imp_mass_amu option.
	*/
	void Options::set_imp_mass_amu(double imp_mass_amu) 
		{m_imp_mass_amu = imp_mass_amu;}
	
	/**
	* Sets the imp_init_charge option.
	*
	* Charge state the impurities start with.
	*
	* @param imp_init_charge New value of the imp_init_charge option.
	*/
	void Options::set_imp_init_charge(int imp_init_charge) 
		{m_imp_init_charge = imp_init_charge;}
	
	/**
	* Sets the imp_num option.
	*
	* Total number of impurity particles to follow.
	*
	* @param imp_num New value of the imp_num option.
	*/
	void Options::set_imp_num(int imp_num) 
		{m_imp_num = imp_num;}
	
	/**
	* Sets the imp_tstart_opt option.
	*
	* How each impurity's starting time is chosen: a single value
	* (imp_tstart_val), uniformly between imp_trange_min and imp_trange_max,
	* or uniformly over the full background time range.
	*
	* Valid values are "single_value", "range" and "full_range". An invalid
	* value prints an error and leaves the option unchanged. Also sets the
	* control integer returned by imp_tstart_opt_int() to 0 for
	* "single_value", 1 for "range", 2 for "full_range". An invalid value
	* leaves the control integer unchanged.
	*
	* @param imp_tstart_opt New value of the imp_tstart_opt option.
	*/
	void Options::set_imp_tstart_opt(std::string imp_tstart_opt) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("imp_tstart_opt", imp_tstart_opt, 
				{"single_value", "range", "full_range"}))
			{
				m_imp_tstart_opt = imp_tstart_opt;
			}

			// Assign control integers
			if (imp_tstart_opt == "single_value") m_imp_tstart_opt_int = 0;
			else if (imp_tstart_opt == "range") m_imp_tstart_opt_int = 1;
			else if (imp_tstart_opt == "full_range") m_imp_tstart_opt_int = 2;
		}

	/**
	* Sets the imp_tstart_val option.
	*
	* Starting time of every impurity when imp_tstart_opt is "single_value"
	* (s).
	*
	* @param imp_tstart_val New value of the imp_tstart_val option.
	*/
	void Options::set_imp_tstart_val(double imp_tstart_val) 
		{m_imp_tstart_val = imp_tstart_val;}
		
	/**
	* Sets the imp_trange_min option.
	*
	* Lower limit of the starting time range when imp_tstart_opt is "range"
	* (s).
	*
	* @param imp_trange_min New value of the imp_trange_min option.
	*/
	void Options::set_imp_trange_min(double imp_trange_min) 
		{m_imp_trange_min = imp_trange_min;}
	
	/**
	* Sets the imp_trange_max option.
	*
	* Upper limit of the starting time range when imp_tstart_opt is "range"
	* (s).
	*
	* @param imp_trange_max New value of the imp_trange_max option.
	*/
	void Options::set_imp_trange_max(double imp_trange_max) 
		{m_imp_trange_max = imp_trange_max;}

	/**
	* Sets the imp_xstart_opt option.
	*
	* How each impurity's starting x coordinate is chosen: a single value
	* (imp_xstart_val), uniformly between imp_xrange_min and imp_xrange_max,
	* or uniformly over the full x range of the background.
	*
	* Valid values are "single_value", "range" and "full_range". An invalid
	* value prints an error and leaves the option unchanged. Also sets the
	* control integer returned by imp_xstart_opt_int() to 0 for
	* "single_value", 1 for "range", 2 for "full_range". An invalid value
	* leaves the control integer unchanged.
	*
	* @param imp_xstart_opt New value of the imp_xstart_opt option.
	*/
	void Options::set_imp_xstart_opt(std::string imp_xstart_opt) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("imp_xstart_opt", imp_xstart_opt, 
				{"single_value", "range", "full_range"}))
			{
				m_imp_xstart_opt = imp_xstart_opt;
			}

			// Assign control integers
			if (imp_xstart_opt == "single_value") m_imp_xstart_opt_int = 0;
			else if (imp_xstart_opt == "range") m_imp_xstart_opt_int = 1;
			else if (imp_xstart_opt == "full_range") m_imp_xstart_opt_int = 2;
		}

	/**
	* Sets the imp_xstart_val option.
	*
	* Starting x coordinate of every impurity when imp_xstart_opt is
	* "single_value".
	*
	* @param imp_xstart_val New value of the imp_xstart_val option.
	*/
	void Options::set_imp_xstart_val(double imp_xstart_val) 
		{m_imp_xstart_val = imp_xstart_val;}
		
	/**
	* Sets the imp_xrange_min option.
	*
	* Lower limit of the starting x range when imp_xstart_opt is "range".
	*
	* @param imp_xrange_min New value of the imp_xrange_min option.
	*/
	void Options::set_imp_xrange_min(double imp_xrange_min) 
		{m_imp_xrange_min = imp_xrange_min;}
	
	/**
	* Sets the imp_xrange_max option.
	*
	* Upper limit of the starting x range when imp_xstart_opt is "range".
	*
	* @param imp_xrange_max New value of the imp_xrange_max option.
	*/
	void Options::set_imp_xrange_max(double imp_xrange_max) 
		{m_imp_xrange_max = imp_xrange_max;}
	
	/**
	* Sets the imp_ystart_opt option.
	*
	* How each impurity's starting y coordinate is chosen: a single value
	* (imp_ystart_val), uniformly between imp_yrange_min and imp_yrange_max,
	* or uniformly over the full y range of the background.
	*
	* Valid values are "single_value", "range" and "full_range". An invalid
	* value prints an error and leaves the option unchanged. Also sets the
	* control integer returned by imp_ystart_opt_int() to 0 for
	* "single_value", 1 for "range", 2 for "full_range". An invalid value
	* leaves the control integer unchanged.
	*
	* @param imp_ystart_opt New value of the imp_ystart_opt option.
	*/
	void Options::set_imp_ystart_opt(std::string imp_ystart_opt) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("imp_ystart_opt", imp_ystart_opt, 
				{"single_value", "range", "full_range"}))
			{
				m_imp_ystart_opt = imp_ystart_opt;
			}

			// Assign control integers
			if (imp_ystart_opt == "single_value") m_imp_ystart_opt_int = 0;
			else if (imp_ystart_opt == "range") m_imp_ystart_opt_int = 1;
			else if (imp_ystart_opt == "full_range") m_imp_ystart_opt_int = 2;
		}
	
	/**
	* Sets the imp_ystart_val option.
	*
	* Starting y coordinate of every impurity when imp_ystart_opt is
	* "single_value".
	*
	* @param imp_ystart_val New value of the imp_ystart_val option.
	*/
	void Options::set_imp_ystart_val(double imp_ystart_val) 
		{m_imp_ystart_val = imp_ystart_val;}

	/**
	* Sets the imp_yrange_min option.
	*
	* Lower limit of the starting y range when imp_ystart_opt is "range".
	*
	* @param imp_yrange_min New value of the imp_yrange_min option.
	*/
	void Options::set_imp_yrange_min(double imp_yrange_min) 
		{m_imp_yrange_min = imp_yrange_min;}
	
	/**
	* Sets the imp_yrange_max option.
	*
	* Upper limit of the starting y range when imp_ystart_opt is "range".
	*
	* @param imp_yrange_max New value of the imp_yrange_max option.
	*/
	void Options::set_imp_yrange_max(double imp_yrange_max) 
		{m_imp_yrange_max = imp_yrange_max;}
	
	/**
	* Sets the imp_zstart_opt option.
	*
	* How each impurity's starting z coordinate is chosen: a single value
	* (imp_zstart_val), uniformly between imp_zrange_min and imp_zrange_max,
	* or uniformly over the full z range of the background.
	*
	* Valid values are "single_value", "range" and "full_range". An invalid
	* value prints an error and leaves the option unchanged. Also sets the
	* control integer returned by imp_zstart_opt_int() to 0 for
	* "single_value", 1 for "range", 2 for "full_range". An invalid value
	* leaves the control integer unchanged.
	*
	* @param imp_zstart_opt New value of the imp_zstart_opt option.
	*/
	void Options::set_imp_zstart_opt(std::string imp_zstart_opt) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("imp_zstart_opt", imp_zstart_opt, 
				{"single_value", "range", "full_range"}))
			{
				m_imp_zstart_opt = imp_zstart_opt;
			}

			// Assign control integers
			if (imp_zstart_opt == "single_value") m_imp_zstart_opt_int = 0;
			else if (imp_zstart_opt == "range") m_imp_zstart_opt_int = 1;
			else if (imp_zstart_opt == "full_range") m_imp_zstart_opt_int = 2;
		}
	
	/**
	* Sets the imp_zstart_val option.
	*
	* Starting z coordinate of every impurity when imp_zstart_opt is
	* "single_value".
	*
	* @param imp_zstart_val New value of the imp_zstart_val option.
	*/
	void Options::set_imp_zstart_val(double imp_zstart_val) 
		{m_imp_zstart_val = imp_zstart_val;}

	/**
	* Sets the imp_zrange_min option.
	*
	* Lower limit of the starting z range when imp_zstart_opt is "range".
	*
	* @param imp_zrange_min New value of the imp_zrange_min option.
	*/
	void Options::set_imp_zrange_min(double imp_zrange_min) 
		{m_imp_zrange_min = imp_zrange_min;}
	
	/**
	* Sets the imp_zrange_max option.
	*
	* Upper limit of the starting z range when imp_zstart_opt is "range".
	*
	* @param imp_zrange_max New value of the imp_zrange_max option.
	*/
	void Options::set_imp_zrange_max(double imp_zrange_max) 
		{m_imp_zrange_max = imp_zrange_max;}
	
	/**
	* Sets the imp_temp_start_opt option.
	*
	* How each impurity's starting temperature is chosen: a single value
	* (imp_temp_start_val) or the local main ion temperature. Starting
	* velocities are sampled from an isotropic Maxwellian at this
	* temperature.
	*
	* Valid values are "single_value" and "main_ion". An invalid value
	* prints an error and leaves the option unchanged. Also sets the control
	* integer returned by imp_temp_start_opt_int() to 0 for "single_value",
	* 1 for "main_ion". An invalid value leaves the control integer
	* unchanged.
	*
	* @param imp_temp_start_opt New value of the imp_temp_start_opt option.
	*/
	void Options::set_imp_temp_start_opt(std::string imp_temp_start_opt) 
	{
		// Only assign if valid option, leaving as default if not
		if (check_input<std::string>("imp_temp_start_opt", 
			imp_temp_start_opt, {"single_value", "main_ion"}))
		{
			m_imp_temp_start_opt = imp_temp_start_opt;
		}

		// Assign control integers
		if (imp_temp_start_opt == "single_value") 
			m_imp_temp_start_opt_int = 0;
		else if (imp_temp_start_opt == "main_ion") 
			m_imp_temp_start_opt_int = 1;
	}
	
	/**
	* Sets the imp_temp_start_val option.
	*
	* Starting temperature of every impurity when imp_temp_start_opt is
	* "single_value" (eV).
	*
	* @param imp_temp_start_val New value of the imp_temp_start_val option.
	*/
	void Options::set_imp_temp_start_val(double imp_temp_start_val) 
		{m_imp_temp_start_val = imp_temp_start_val;}

	/**
	* Sets the imp_collisions option.
	*
	* Collision model for the impurities, either off or the Nanbu model.
	*
	* Valid values are "off" and "nanbu". An invalid value prints an error
	* and leaves the option unchanged. Also sets the control integer
	* returned by imp_collisions_int() to 0 for "off", 1 for "nanbu". An
	* invalid value leaves the control integer unchanged.
	*
	* @param imp_collisions New value of the imp_collisions option.
	*/
	void Options::set_imp_collisions(std::string imp_collisions) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("imp_collision", imp_collisions, 
				{"off", "nanbu"}))
			{
				m_imp_collisions = imp_collisions;
			}

			// Assign control integers
			if (imp_collisions == "off") m_imp_collisions_int = 0;
			else if (imp_collisions == "nanbu") m_imp_collisions_int = 1;
		}

	/**
	* Sets the var_red_split option.
	*
	* When particles are split for variance reduction: never, at
	* ionization/recombination events, or at collisions.
	*
	* Valid values are "off", "iz_rec" and "coll". An invalid value prints
	* an error and leaves the option unchanged. Also sets the control
	* integer returned by var_red_split_int() to 0 for "off", 1 for
	* "iz_rec", 2 for "coll". An invalid value leaves the control integer
	* unchanged.
	*
	* @param var_red_split New value of the var_red_split option.
	*/
	void Options::set_var_red_split(std::string var_red_split) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("var_red_split", var_red_split, 
				{"off", "iz_rec", "coll"}))
			{
				m_var_red_split = var_red_split;
			}

			// Assign control integers
			if (var_red_split == "off") m_var_red_split_int = 0;
			else if (var_red_split == "iz_rec") m_var_red_split_int = 1;
			else if (var_red_split == "coll") m_var_red_split_int = 2;
		}

	/**
	* Sets the var_red_import option.
	*
	* Method used to calculate particle importance for variance reduction.
	*
	* Valid values are "median", "exp_dist" and "exp_time". An invalid value
	* prints an error and leaves the option unchanged. Also sets the control
	* integer returned by var_red_import_int() to 0 for "median", 1 for
	* "exp_dist", 2 for "exp_time". An invalid value leaves the control
	* integer unchanged.
	*
	* @param var_red_import New value of the var_red_import option.
	*/
	void Options::set_var_red_import(std::string var_red_import) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("var_red_import", var_red_import, 
				{"median", "exp_dist", "exp_time"}))
			{
				m_var_red_import = var_red_import;
			}

			// Assign control integers
			if (var_red_import == "median") m_var_red_import_int = 0;
			else if (var_red_import == "exp_dist") m_var_red_import_int = 1;
			else if (var_red_import == "exp_time") m_var_red_import_int = 2;
		}

	/**
	* Sets the var_red_freq option.
	*
	* How often variance reduction is applied.
	*
	* @param var_red_freq New value of the var_red_freq option.
	*/
	void Options::set_var_red_freq(double var_red_freq) 
		{m_var_red_freq = var_red_freq;}
	
	/**
	* Sets the var_red_min_weight option.
	*
	* Minimum particle weight for variance reduction.
	*
	* @param var_red_min_weight New value of the var_red_min_weight option.
	*/
	void Options::set_var_red_min_weight(double var_red_min_weight) 
		{m_var_red_min_weight = var_red_min_weight;}
	
	/**
	* Sets the var_red_med_mod option.
	*
	* Modifier applied to the median importance for variance reduction.
	*
	* @param var_red_med_mod New value of the var_red_med_mod option.
	*/
	void Options::set_var_red_med_mod(double var_red_med_mod) 
		{m_var_red_med_mod = var_red_med_mod;}
	
	/**
	* Sets the var_red_rusrol option.
	*
	* Whether to use Russian roulette for variance reduction.
	*
	* Valid values are "off" and "on". An invalid value prints an error and
	* leaves the option unchanged. Also sets the control integer returned by
	* var_red_rusrol_int() to 0 for "off", 1 for "on". An invalid value
	* leaves the control integer unchanged.
	*
	* @param var_red_rusrol New value of the var_red_rusrol option.
	*/
	void Options::set_var_red_rusrol(std::string var_red_rusrol) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("var_red_rusrol", var_red_rusrol, 
				{"off", "on"}))
			{
				m_var_red_rusrol = var_red_rusrol;
			}

			// Assign control integers
			if (var_red_rusrol == "off") m_var_red_rusrol_int = 0;
			else if (var_red_rusrol == "on") m_var_red_rusrol_int = 1;
		}

	/**
	* Sets the var_red_rusrol_prob option.
	*
	* Probability used for Russian roulette.
	*
	* @param var_red_rusrol_prob New value of the var_red_rusrol_prob
	*                            option.
	*/
	void Options::set_var_red_rusrol_prob(double var_red_rusrol_prob) 
		{m_var_red_rusrol_prob = var_red_rusrol_prob;}

	/**
	* Sets the imp_time_step_opt option.
	*
	* Whether the impurity time step is constant or variable.
	*
	* Valid values are "constant" and "variable". An invalid value prints an
	* error and leaves the option unchanged. Also sets the control integer
	* returned by imp_time_step_opt_int() to 0 for "constant", 1 for
	* "variable". An invalid value leaves the control integer unchanged.
	*
	* @param imp_time_step_opt New value of the imp_time_step_opt option.
	*/
	void Options::set_imp_time_step_opt(std::string imp_time_step_opt) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("imp_time_step_opt", imp_time_step_opt, 
				{"constant", "variable"}))
			{
				m_imp_time_step_opt = imp_time_step_opt;
			}

			// Assign control integers
			if (imp_time_step_opt == "constant") m_imp_time_step_opt_int = 0;
			else if (imp_time_step_opt == "variable") 
				m_imp_time_step_opt_int = 1;
		}
	
	/**
	* Sets the imp_time_step option.
	*
	* Impurity time step (s).
	*
	* @param imp_time_step New value of the imp_time_step option.
	*/
	void Options::set_imp_time_step(double imp_time_step) 
		{m_imp_time_step = imp_time_step;}
	
	/**
	* Sets the imp_time_step_min option.
	*
	* Minimum impurity time step when the time step is variable (s).
	*
	* @param imp_time_step_min New value of the imp_time_step_min option.
	*/
	void Options::set_imp_time_step_min(double imp_time_step_min) 
		{m_imp_time_step_min = imp_time_step_min;}
	
	/**
	* Sets the imp_source_scale_fact option.
	*
	* Impurity source scaling factor (particles/s), i.e., how many real
	* particles a single simulated particle represents. Used to convert the
	* Monte Carlo counts into a density.
	*
	* @param imp_source_scale_fact New value of the imp_source_scale_fact
	*                              option.
	*/
	void Options::set_imp_source_scale_fact(double imp_source_scale_fact) 
		{m_imp_source_scale_fact = imp_source_scale_fact;}
	
	/**
	* Sets the imp_xbound_buffer option.
	*
	* Buffer distance at the x boundaries. With an absorbing boundary, a
	* particle is absorbed once it is within this distance of the boundary.
	*
	* @param imp_xbound_buffer New value of the imp_xbound_buffer option.
	*/
	void Options::set_imp_xbound_buffer(double imp_xbound_buffer) 
		{m_imp_xbound_buffer = imp_xbound_buffer;}

	/**
	* Sets the imp_ybound_buffer option.
	*
	* Buffer distance at the y boundaries. With an absorbing boundary, a
	* particle is absorbed once it is within this distance of the boundary.
	*
	* @param imp_ybound_buffer New value of the imp_ybound_buffer option.
	*/
	void Options::set_imp_ybound_buffer(double imp_ybound_buffer) 
		{m_imp_ybound_buffer = imp_ybound_buffer;}

	/**
	* Sets the imp_zbound_buffer option.
	*
	* Buffer distance at the z boundaries. With an absorbing boundary, a
	* particle is absorbed once it is within this distance of the boundary.
	*
	* @param imp_zbound_buffer New value of the imp_zbound_buffer option.
	*/
	void Options::set_imp_zbound_buffer(double imp_zbound_buffer) 
		{m_imp_zbound_buffer = imp_zbound_buffer;}
	
	/**
	* Sets the imp_iz_recomb option.
	*
	* Whether impurity ionization and recombination are on.
	*
	* Valid values are "off" and "on". An invalid value prints an error and
	* leaves the option unchanged. Also sets the control integer returned by
	* imp_iz_recomb_int() to 0 for "off", 1 for "on". An invalid value
	* leaves the control integer unchanged.
	*
	* @param imp_iz_recomb New value of the imp_iz_recomb option.
	*/
	void Options::set_imp_iz_recomb(std::string imp_iz_recomb) 
		{
			// Only assign if valid option, leaving as default if not
			if (check_input<std::string>("imp_iz_recomb", imp_iz_recomb, 
				{"off", "on"}))
			{
				m_imp_iz_recomb = imp_iz_recomb;
			}

			// Set control integers
			if (imp_iz_recomb == "off") m_imp_iz_recomb_int = 0;
			else if (imp_iz_recomb == "on") m_imp_iz_recomb_int = 1;
		}

	/**
	* Sets the imp_max_time option.
	*
	* Maximum time a particle is followed for (s). A particle that exceeds
	* it is killed and counted in a warning at the end of the run.
	*
	* @param imp_max_time New value of the imp_max_time option.
	*/
	void Options::set_imp_max_time(double imp_max_time) 
		{m_imp_max_time = imp_max_time;}
	
	/**
	* Sets the print_interval option.
	*
	* Number of progress updates printed while particles are followed. For
	* example, 10 prints an update every 10% of particles.
	*
	* @param print_interval New value of the print_interval option.
	*/
	void Options::set_print_interval(int print_interval) 
		{m_print_interval = print_interval;}
	
	/**
	* Sets the openadas_root option.
	*
	* Directory containing the OpenADAS rate coefficient files.
	*
	* @param openadas_root New value of the openadas_root option.
	*/
	void Options::set_openadas_root(std::string openadas_root) 
		{m_openadas_root = openadas_root;}
	
	/**
	* Sets the openadas_year option.
	*
	* Year of the OpenADAS data files to use.
	*
	* @param openadas_year New value of the openadas_year option.
	*/
	void Options::set_openadas_year(int openadas_year) 
		{m_openadas_year = openadas_year;}
	
	// Accessor definitions

	/**
	* Returns the case_name option.
	*
	* Name of the simulation case. The results are saved to <case_name>.nc.
	*
	* @return Value of the case_name option.
	*/
	const std::string& Options::case_name() const 
		{return m_case_name;}

	/**
	* Returns the use_gpu option.
	*
	* Whether to run the simulation on GPUs.
	*
	* @return Value of the use_gpu option.
	*/
	const std::string& Options::use_gpu() const 
		{return m_use_gpu;}

	/**
	* Returns the seed option.
	*
	* Base seed for the random number generators. A negative value seeds
	* them from the clock so each run is different. The MPI rank is added to
	* the seed so each rank gets a different sequence.
	*
	* @return Value of the seed option.
	*/
	int Options::seed() const 
		{return m_seed;}

	/**
	* Returns the slot_cap option.
	*
	* Maximum number of particle slots, i.e., the most particles followed at
	* once. The number of slots used is the smaller of this and imp_num.
	*
	* @return Value of the slot_cap option.
	*/
	int Options::slot_cap() const 
		{return m_slot_cap;}

	/**
	* Returns the bkg_source option.
	*
	* Source of the background plasma, either a built-in test background or
	* a Gkeyll simulation.
	*
	* @return Value of the bkg_source option.
	*/
	const std::string& Options::bkg_source() const 
		{return m_bkg_source;}

	/**
	* Returns the test_opt option.
	*
	* Which built-in test case to run when using the test background.
	*
	* @return Value of the test_opt option.
	*/
	const std::string& Options::test_opt() const 
		{return m_test_opt;}

	/**
	* Returns the save_track option.
	*
	* Whether to save particle tracks to the output file.
	*
	* @return Value of the save_track option.
	*/
	const std::string& Options::save_track() const 
		{return m_save_track;}

	/**
	* Returns the gkyl_dir option.
	*
	* Directory containing the Gkeyll simulation files.
	*
	* @return Value of the gkyl_dir option.
	*/
	const std::string& Options::gkyl_dir() const 
		{return m_gkyl_dir;}

	/**
	* Returns the gkyl_casename option.
	*
	* Name of the Gkeyll case, i.e., the prefix of its output files.
	*
	* @return Value of the gkyl_casename option.
	*/
	const std::string& Options::gkyl_casename() const 
		{return m_gkyl_casename;}

	/**
	* Returns the gkyl_frame_start option.
	*
	* First Gkeyll frame to load.
	*
	* @return Value of the gkyl_frame_start option.
	*/
	int Options::gkyl_frame_start() const 
		{return m_gkyl_frame_start;}

	/**
	* Returns the gkyl_frame_end option.
	*
	* Last Gkeyll frame to load. The range from gkyl_frame_start to
	* gkyl_frame_end is inclusive.
	*
	* @return Value of the gkyl_frame_end option.
	*/
	int Options::gkyl_frame_end() const 
		{return m_gkyl_frame_end;}

	/**
	* Returns the gkyl_elec_name option.
	*
	* Name of the electron species in the Gkeyll simulation.
	*
	* @return Value of the gkyl_elec_name option.
	*/
	const std::string& Options::gkyl_elec_name() const 
		{return m_gkyl_elec_name;}

	/**
	* Returns the gkyl_ion_name option.
	*
	* Name of the main ion species in the Gkeyll simulation.
	*
	* @return Value of the gkyl_ion_name option.
	*/
	const std::string& Options::gkyl_ion_name() const 
		{return m_gkyl_ion_name;}

	/**
	* Returns the gkyl_elec_mass_amu option.
	*
	* Electron mass used in the Gkeyll simulation (amu).
	*
	* @return Value of the gkyl_elec_mass_amu option.
	*/
	double Options::gkyl_elec_mass_amu() const 
		{return m_gkyl_elec_mass_amu;}

	/**
	* Returns the gkyl_ion_mass_amu option.
	*
	* Main ion mass used in the Gkeyll simulation (amu).
	*
	* @return Value of the gkyl_ion_mass_amu option.
	*/
	double Options::gkyl_ion_mass_amu() const 
		{return m_gkyl_ion_mass_amu;}

	/**
	* Returns the gkyl_file_type option.
	*
	* File format of the Gkeyll output.
	*
	* @return Value of the gkyl_file_type option.
	*/
	const std::string& Options::gkyl_file_type() const 
		{return m_gkyl_file_type;}

	/**
	* Returns the gkyl_moment_type option.
	*
	* Type of moment files in the Gkeyll output. Collisions need
	* "bimaxwellian", which provides the ion parallel flow and perpendicular
	* velocity used to calculate the background flow.
	*
	* @return Value of the gkyl_moment_type option.
	*/
	const std::string& Options::gkyl_moment_type() const 
		{return m_gkyl_moment_type;}

	/**
	* Returns the min_ne option.
	*
	* Lower limit on the electron density (m^-3).
	*
	* @return Value of the min_ne option.
	*/
	double Options::min_ne() const 
		{return m_min_ne;}

	/**
	* Returns the min_te option.
	*
	* Lower limit on the electron temperature (eV).
	*
	* @return Value of the min_te option.
	*/
	double Options::min_te() const 
		{return m_min_te;}

	/**
	* Returns the min_ti option.
	*
	* Lower limit on the ion temperature (eV), applied when calculating
	* collisions with the background ions.
	*
	* @return Value of the min_ti option.
	*/
	double Options::min_ti() const 
		{return m_min_ti;}

	/**
	* Returns the min_xbound_type option.
	*
	* Boundary condition at the minimum x boundary.
	*
	* @return Value of the min_xbound_type option.
	*/
	const std::string& Options::min_xbound_type() const 
		{return m_min_xbound_type;}

	/**
	* Returns the max_xbound_type option.
	*
	* Boundary condition at the maximum x boundary.
	*
	* @return Value of the max_xbound_type option.
	*/
	const std::string& Options::max_xbound_type() const 
		{return m_max_xbound_type;}

	/**
	* Returns the min_ybound_type option.
	*
	* Boundary condition at the minimum y boundary.
	*
	* @return Value of the min_ybound_type option.
	*/
	const std::string& Options::min_ybound_type() const 
		{return m_min_ybound_type;}

	/**
	* Returns the max_ybound_type option.
	*
	* Boundary condition at the maximum y boundary.
	*
	* @return Value of the max_ybound_type option.
	*/
	const std::string& Options::max_ybound_type() const 
		{return m_max_ybound_type;}

	/**
	* Returns the min_zbound_type option.
	*
	* Boundary condition at the minimum z boundary.
	*
	* @return Value of the min_zbound_type option.
	*/
	const std::string& Options::min_zbound_type() const 
		{return m_min_zbound_type;}

	/**
	* Returns the max_zbound_type option.
	*
	* Boundary condition at the maximum z boundary.
	*
	* @return Value of the max_zbound_type option.
	*/
	const std::string& Options::max_zbound_type() const 
		{return m_max_zbound_type;}

	/**
	* Returns the lcfs_x option.
	*
	* x coordinate of the last closed flux surface.
	*
	* @return Value of the lcfs_x option.
	*/
	double Options::lcfs_x() const 
		{return m_lcfs_x;}

	/**
	* Returns the lcfs_y option.
	*
	* y coordinate of the last closed flux surface.
	*
	* @return Value of the lcfs_y option.
	*/
	double Options::lcfs_y() const 
		{return m_lcfs_y;}

	/**
	* Returns the lcfs_z option.
	*
	* z coordinate of the last closed flux surface.
	*
	* @return Value of the lcfs_z option.
	*/
	double Options::lcfs_z() const 
		{return m_lcfs_z;}

	/**
	* Returns the sep_x_bc_xp_z1 option.
	*
	* First z value for the separatrix x boundary condition near the
	* X-point.
	*
	* @return Value of the sep_x_bc_xp_z1 option.
	*/
	double Options::sep_x_bc_xp_z1() const 
		{return m_sep_x_bc_xp_z1;}

	/**
	* Returns the sep_x_bc_xp_z2 option.
	*
	* Second z value for the separatrix x boundary condition near the
	* X-point.
	*
	* @return Value of the sep_x_bc_xp_z2 option.
	*/
	double Options::sep_x_bc_xp_z2() const 
		{return m_sep_x_bc_xp_z2;}

	/**
	* Returns the imp_xbound_buffer option.
	*
	* Buffer distance at the x boundaries. With an absorbing boundary, a
	* particle is absorbed once it is within this distance of the boundary.
	*
	* @return Value of the imp_xbound_buffer option.
	*/
	double Options::imp_xbound_buffer() const 
		{return m_imp_xbound_buffer;}

	/**
	* Returns the imp_ybound_buffer option.
	*
	* Buffer distance at the y boundaries. With an absorbing boundary, a
	* particle is absorbed once it is within this distance of the boundary.
	*
	* @return Value of the imp_ybound_buffer option.
	*/
	double Options::imp_ybound_buffer() const 
		{return m_imp_ybound_buffer;}

	/**
	* Returns the imp_zbound_buffer option.
	*
	* Buffer distance at the z boundaries. With an absorbing boundary, a
	* particle is absorbed once it is within this distance of the boundary.
	*
	* @return Value of the imp_zbound_buffer option.
	*/
	double Options::imp_zbound_buffer() const 
		{return m_imp_zbound_buffer;}

	/**
	* Returns the imp_atom_num option.
	*
	* Atomic number of the impurity.
	*
	* @return Value of the imp_atom_num option.
	*/
	int Options::imp_atom_num() const 
		{return m_imp_atom_num;}

	/**
	* Returns the imp_mass_amu option.
	*
	* Mass of the impurity (amu).
	*
	* @return Value of the imp_mass_amu option.
	*/
	double Options::imp_mass_amu() const 
		{return m_imp_mass_amu;}

	/**
	* Returns the imp_init_charge option.
	*
	* Charge state the impurities start with.
	*
	* @return Value of the imp_init_charge option.
	*/
	int Options::imp_init_charge() const 
		{return m_imp_init_charge;}

	/**
	* Returns the imp_num option.
	*
	* Total number of impurity particles to follow.
	*
	* @return Value of the imp_num option.
	*/
	int Options::imp_num() const
		{return m_imp_num;}

	/**
	* Returns the imp_tstart_opt option.
	*
	* How each impurity's starting time is chosen: a single value
	* (imp_tstart_val), uniformly between imp_trange_min and imp_trange_max,
	* or uniformly over the full background time range.
	*
	* @return Value of the imp_tstart_opt option.
	*/
	const std::string& Options::imp_tstart_opt() const 
		{return m_imp_tstart_opt;}

	/**
	* Returns the imp_tstart_val option.
	*
	* Starting time of every impurity when imp_tstart_opt is "single_value"
	* (s).
	*
	* @return Value of the imp_tstart_val option.
	*/
	double Options::imp_tstart_val() const 
		{return m_imp_tstart_val;}

	/**
	* Returns the imp_trange_min option.
	*
	* Lower limit of the starting time range when imp_tstart_opt is "range"
	* (s).
	*
	* @return Value of the imp_trange_min option.
	*/
	double Options::imp_trange_min() const 
		{return m_imp_trange_min;}

	/**
	* Returns the imp_trange_max option.
	*
	* Upper limit of the starting time range when imp_tstart_opt is "range"
	* (s).
	*
	* @return Value of the imp_trange_max option.
	*/
	double Options::imp_trange_max() const 
		{return m_imp_trange_max;}

	/**
	* Returns the imp_xstart_opt option.
	*
	* How each impurity's starting x coordinate is chosen: a single value
	* (imp_xstart_val), uniformly between imp_xrange_min and imp_xrange_max,
	* or uniformly over the full x range of the background.
	*
	* @return Value of the imp_xstart_opt option.
	*/
	const std::string& Options::imp_xstart_opt() const 
		{return m_imp_xstart_opt;}

	/**
	* Returns the imp_xstart_val option.
	*
	* Starting x coordinate of every impurity when imp_xstart_opt is
	* "single_value".
	*
	* @return Value of the imp_xstart_val option.
	*/
	double Options::imp_xstart_val() const 
		{return m_imp_xstart_val;}

	/**
	* Returns the imp_xrange_min option.
	*
	* Lower limit of the starting x range when imp_xstart_opt is "range".
	*
	* @return Value of the imp_xrange_min option.
	*/
	double Options::imp_xrange_min() const 
		{return m_imp_xrange_min;}

	/**
	* Returns the imp_xrange_max option.
	*
	* Upper limit of the starting x range when imp_xstart_opt is "range".
	*
	* @return Value of the imp_xrange_max option.
	*/
	double Options::imp_xrange_max() const 
		{return m_imp_xrange_max;}

	/**
	* Returns the imp_ystart_opt option.
	*
	* How each impurity's starting y coordinate is chosen: a single value
	* (imp_ystart_val), uniformly between imp_yrange_min and imp_yrange_max,
	* or uniformly over the full y range of the background.
	*
	* @return Value of the imp_ystart_opt option.
	*/
	const std::string& Options::imp_ystart_opt() const 
		{return m_imp_ystart_opt;}

	/**
	* Returns the imp_ystart_val option.
	*
	* Starting y coordinate of every impurity when imp_ystart_opt is
	* "single_value".
	*
	* @return Value of the imp_ystart_val option.
	*/
	double Options::imp_ystart_val() const 
		{return m_imp_ystart_val;}

	/**
	* Returns the imp_yrange_min option.
	*
	* Lower limit of the starting y range when imp_ystart_opt is "range".
	*
	* @return Value of the imp_yrange_min option.
	*/
	double Options::imp_yrange_min() const 
		{return m_imp_yrange_min;}

	/**
	* Returns the imp_yrange_max option.
	*
	* Upper limit of the starting y range when imp_ystart_opt is "range".
	*
	* @return Value of the imp_yrange_max option.
	*/
	double Options::imp_yrange_max() const 
		{return m_imp_yrange_max;}

	/**
	* Returns the imp_zstart_opt option.
	*
	* How each impurity's starting z coordinate is chosen: a single value
	* (imp_zstart_val), uniformly between imp_zrange_min and imp_zrange_max,
	* or uniformly over the full z range of the background.
	*
	* @return Value of the imp_zstart_opt option.
	*/
	const std::string& Options::imp_zstart_opt() const 
		{return m_imp_zstart_opt;}

	/**
	* Returns the imp_zstart_val option.
	*
	* Starting z coordinate of every impurity when imp_zstart_opt is
	* "single_value".
	*
	* @return Value of the imp_zstart_val option.
	*/
	double Options::imp_zstart_val() const 
		{return m_imp_zstart_val;}

	/**
	* Returns the imp_zrange_min option.
	*
	* Lower limit of the starting z range when imp_zstart_opt is "range".
	*
	* @return Value of the imp_zrange_min option.
	*/
	double Options::imp_zrange_min() const 
		{return m_imp_zrange_min;}

	/**
	* Returns the imp_zrange_max option.
	*
	* Upper limit of the starting z range when imp_zstart_opt is "range".
	*
	* @return Value of the imp_zrange_max option.
	*/
	double Options::imp_zrange_max() const 
		{return m_imp_zrange_max;}

	/**
	* Returns the imp_temp_start_opt option.
	*
	* How each impurity's starting temperature is chosen: a single value
	* (imp_temp_start_val) or the local main ion temperature. Starting
	* velocities are sampled from an isotropic Maxwellian at this
	* temperature.
	*
	* @return Value of the imp_temp_start_opt option.
	*/
	const std::string& Options::imp_temp_start_opt() const 
		{return m_imp_temp_start_opt;}

	/**
	* Returns the imp_temp_start_val option.
	*
	* Starting temperature of every impurity when imp_temp_start_opt is
	* "single_value" (eV).
	*
	* @return Value of the imp_temp_start_val option.
	*/
	double Options::imp_temp_start_val() const 
		{return m_imp_temp_start_val;}

	/**
	* Returns the imp_collisions option.
	*
	* Collision model for the impurities, either off or the Nanbu model.
	*
	* @return Value of the imp_collisions option.
	*/
	const std::string& Options::imp_collisions() const 
		{return m_imp_collisions;}

	/**
	* Returns the var_red_split option.
	*
	* When particles are split for variance reduction: never, at
	* ionization/recombination events, or at collisions.
	*
	* @return Value of the var_red_split option.
	*/
	const std::string& Options::var_red_split() const 
		{return m_var_red_split;}

	/**
	* Returns the var_red_import option.
	*
	* Method used to calculate particle importance for variance reduction.
	*
	* @return Value of the var_red_import option.
	*/
	const std::string& Options::var_red_import() const 
		{return m_var_red_import;}

	/**
	* Returns the var_red_freq option.
	*
	* How often variance reduction is applied.
	*
	* @return Value of the var_red_freq option.
	*/
	double Options::var_red_freq() const 
		{return m_var_red_freq;}

	/**
	* Returns the var_red_min_weight option.
	*
	* Minimum particle weight for variance reduction.
	*
	* @return Value of the var_red_min_weight option.
	*/
	double Options::var_red_min_weight() const 
		{return m_var_red_min_weight;}

	/**
	* Returns the var_red_med_mod option.
	*
	* Modifier applied to the median importance for variance reduction.
	*
	* @return Value of the var_red_med_mod option.
	*/
	double Options::var_red_med_mod() const 
		{return m_var_red_med_mod;}

	/**
	* Returns the var_red_rusrol option.
	*
	* Whether to use Russian roulette for variance reduction.
	*
	* @return Value of the var_red_rusrol option.
	*/
	const std::string& Options::var_red_rusrol() const 
		{return m_var_red_rusrol;}

	/**
	* Returns the var_red_rusrol_prob option.
	*
	* Probability used for Russian roulette.
	*
	* @return Value of the var_red_rusrol_prob option.
	*/
	double Options::var_red_rusrol_prob() const 
		{return m_var_red_rusrol_prob;}

	/**
	* Returns the imp_time_step_opt option.
	*
	* Whether the impurity time step is constant or variable.
	*
	* @return Value of the imp_time_step_opt option.
	*/
	const std::string& Options::imp_time_step_opt() const 
		{return m_imp_time_step_opt;}

	/**
	* Returns the imp_time_step option.
	*
	* Impurity time step (s).
	*
	* @return Value of the imp_time_step option.
	*/
	double Options::imp_time_step() const 
		{return m_imp_time_step;}

	/**
	* Returns the imp_time_step_min option.
	*
	* Minimum impurity time step when the time step is variable (s).
	*
	* @return Value of the imp_time_step_min option.
	*/
	double Options::imp_time_step_min() const 
		{return m_imp_time_step_min;}

	/**
	* Returns the imp_source_scale_fact option.
	*
	* Impurity source scaling factor (particles/s), i.e., how many real
	* particles a single simulated particle represents. Used to convert the
	* Monte Carlo counts into a density.
	*
	* @return Value of the imp_source_scale_fact option.
	*/
	double Options::imp_source_scale_fact() const 
		{return m_imp_source_scale_fact;}

	/**
	* Returns the imp_iz_recomb option.
	*
	* Whether impurity ionization and recombination are on.
	*
	* @return Value of the imp_iz_recomb option.
	*/
	const std::string& Options::imp_iz_recomb() const 
		{return m_imp_iz_recomb;}

	/**
	* Returns the imp_max_time option.
	*
	* Maximum time a particle is followed for (s). A particle that exceeds
	* it is killed and counted in a warning at the end of the run.
	*
	* @return Value of the imp_max_time option.
	*/
	double Options::imp_max_time() const 
		{return m_imp_max_time;}

	/**
	* Returns the print_interval option.
	*
	* Number of progress updates printed while particles are followed. For
	* example, 10 prints an update every 10% of particles.
	*
	* @return Value of the print_interval option.
	*/
	int Options::print_interval() const 
		{return m_print_interval;}

	/**
	* Returns the openadas_root option.
	*
	* Directory containing the OpenADAS rate coefficient files.
	*
	* @return Value of the openadas_root option.
	*/
	const std::string& Options::openadas_root() const 
		{return m_openadas_root;}

	/**
	* Returns the openadas_year option.
	*
	* Year of the OpenADAS data files to use.
	*
	* @return Value of the openadas_year option.
	*/
	int Options::openadas_year() const 
		{return m_openadas_year;}
	
	/**
	* Returns a control integer after checking that it was set.
	*
	* Control integers start at -1, meaning not set. If control_int is still
	* negative, an error is printed. This usually means the option's setter
	* is missing from initialize_control_ints. The value is returned either
	* way.
	*
	* @param var_name Name of the option the control integer belongs to (for
	*                 the error message).
	* @param control_int Control integer to check.
	* @return The control integer, unchanged.
	*/
	int Options::get_control_int(const std::string& var_name, 
		const int control_int) const
	{
		if (control_int < 0)
		{
			// If you encounter this error, then it probably means the control
			// integer was not being set at the beginning within 
			// initialize_control_ints. Easy mistake to make.
			std::cerr << "Error! Control integer for " << var_name 
				<< " was never set! This is a programming error, contact a "
				<< "developer. See initialize_control_ints for details.\n";
		}
		return control_int;
	}

	// Accessors for internal control variables

	/**
	* Returns the control integer for the use_gpu option.
	*
	* The control integer is an integer version of use_gpu that is faster to
	* check than the string: 0 = "off", 1 = "cuda". It is fetched through
	* get_control_int, which prints an error if it was never set.
	*
	* @return Control integer for use_gpu, or -1 if it was never set.
	*/
	int Options::use_gpu_int() const
		{return get_control_int("use_gpu", m_use_gpu_int);}

	/**
	* Returns the control integer for the bkg_source option.
	*
	* The control integer is an integer version of bkg_source that is faster
	* to check than the string: 0 = "test", 1 = "gkeyll". It is fetched
	* through get_control_int, which prints an error if it was never set.
	*
	* @return Control integer for bkg_source, or -1 if it was never set.
	*/
	int Options::bkg_source_int() const
		{return get_control_int("bkg_source", m_bkg_source_int);}

	/**
	* Returns the control integer for the test_opt option.
	*
	* The control integer is an integer version of test_opt that is faster
	* to check than the string: 0 = "gyrate", 1 = "exb", 2 = "gradb", 3 =
	* "polarization", 4 = "curvature", 5 = "friction_force". It is fetched
	* through get_control_int, which prints an error if it was never set.
	*
	* @return Control integer for test_opt, or -1 if it was never set.
	*/
	int Options::test_opt_int() const
		{return get_control_int("test_opt", m_test_opt_int);}

	/**
	* Returns the control integer for the save_track option.
	*
	* The control integer is an integer version of save_track that is faster
	* to check than the string: 0 = "off", 1 = "on". It is fetched through
	* get_control_int, which prints an error if it was never set.
	*
	* @return Control integer for save_track, or -1 if it was never set.
	*/
	int Options::save_track_int() const
		{return get_control_int("save_track", m_save_track_int);}

	/**
	* Returns the control integer for the imp_tstart_opt option.
	*
	* The control integer is an integer version of imp_tstart_opt that is
	* faster to check than the string: 0 = "single_value", 1 = "range", 2 =
	* "full_range". It is fetched through get_control_int, which prints an
	* error if it was never set.
	*
	* @return Control integer for imp_tstart_opt, or -1 if it was never set.
	*/
	int Options::imp_tstart_opt_int() const
		{return get_control_int("imp_tstart_opt", m_imp_tstart_opt_int);}

	/**
	* Returns the control integer for the imp_xstart_opt option.
	*
	* The control integer is an integer version of imp_xstart_opt that is
	* faster to check than the string: 0 = "single_value", 1 = "range", 2 =
	* "full_range". It is fetched through get_control_int, which prints an
	* error if it was never set.
	*
	* @return Control integer for imp_xstart_opt, or -1 if it was never set.
	*/
	int Options::imp_xstart_opt_int() const
		{return get_control_int("imp_xstart_opt", m_imp_xstart_opt_int);}

	/**
	* Returns the control integer for the imp_ystart_opt option.
	*
	* The control integer is an integer version of imp_ystart_opt that is
	* faster to check than the string: 0 = "single_value", 1 = "range", 2 =
	* "full_range". It is fetched through get_control_int, which prints an
	* error if it was never set.
	*
	* @return Control integer for imp_ystart_opt, or -1 if it was never set.
	*/
	int Options::imp_ystart_opt_int() const
		{return get_control_int("imp_ystart_opt", m_imp_ystart_opt_int);}

	/**
	* Returns the control integer for the imp_zstart_opt option.
	*
	* The control integer is an integer version of imp_zstart_opt that is
	* faster to check than the string: 0 = "single_value", 1 = "range", 2 =
	* "full_range". It is fetched through get_control_int, which prints an
	* error if it was never set.
	*
	* @return Control integer for imp_zstart_opt, or -1 if it was never set.
	*/
	int Options::imp_zstart_opt_int() const
		{return get_control_int("imp_zstart_opt", m_imp_zstart_opt_int);}

	/**
	* Returns the control integer for the imp_temp_start_opt option.
	*
	* The control integer is an integer version of imp_temp_start_opt that
	* is faster to check than the string: 0 = "single_value", 1 =
	* "main_ion". It is fetched through get_control_int, which prints an
	* error if it was never set.
	*
	* @return Control integer for imp_temp_start_opt, or -1 if it was never
	*         set.
	*/
	int Options::imp_temp_start_opt_int() const
		{return get_control_int("imp_temp_start_opt", m_imp_temp_start_opt_int);}

	/**
	* Returns the control integer for the imp_collisions option.
	*
	* The control integer is an integer version of imp_collisions that is
	* faster to check than the string: 0 = "off", 1 = "nanbu". It is fetched
	* through get_control_int, which prints an error if it was never set.
	*
	* @return Control integer for imp_collisions, or -1 if it was never set.
	*/
	int Options::imp_collisions_int() const
		{return get_control_int("imp_collisions", m_imp_collisions_int);}

	/**
	* Returns the control integer for the var_red_split option.
	*
	* The control integer is an integer version of var_red_split that is
	* faster to check than the string: 0 = "off", 1 = "iz_rec", 2 = "coll".
	* It is fetched through get_control_int, which prints an error if it was
	* never set.
	*
	* @return Control integer for var_red_split, or -1 if it was never set.
	*/
	int Options::var_red_split_int() const
		{return get_control_int("var_red_split", m_var_red_split_int);}

	/**
	* Returns the control integer for the var_red_import option.
	*
	* The control integer is an integer version of var_red_import that is
	* faster to check than the string: 0 = "median", 1 = "exp_dist", 2 =
	* "exp_time". It is fetched through get_control_int, which prints an
	* error if it was never set.
	*
	* @return Control integer for var_red_import, or -1 if it was never set.
	*/
	int Options::var_red_import_int() const
		{return get_control_int("var_red_import", m_var_red_import_int);}

	/**
	* Returns the control integer for the var_red_rusrol option.
	*
	* The control integer is an integer version of var_red_rusrol that is
	* faster to check than the string: 0 = "off", 1 = "on". It is fetched
	* through get_control_int, which prints an error if it was never set.
	*
	* @return Control integer for var_red_rusrol, or -1 if it was never set.
	*/
	int Options::var_red_rusrol_int() const
		{return get_control_int("var_red_rusrol", m_var_red_rusrol_int);}

	/**
	* Returns the control integer for the imp_time_step_opt option.
	*
	* The control integer is an integer version of imp_time_step_opt that is
	* faster to check than the string: 0 = "constant", 1 = "variable". It is
	* fetched through get_control_int, which prints an error if it was never
	* set.
	*
	* @return Control integer for imp_time_step_opt, or -1 if it was never
	*         set.
	*/
	int Options::imp_time_step_opt_int() const
		{return get_control_int("imp_time_step_opt", m_imp_time_step_opt_int);}

	/**
	* Returns the control integer for the imp_iz_recomb option.
	*
	* The control integer is an integer version of imp_iz_recomb that is
	* faster to check than the string: 0 = "off", 1 = "on". It is fetched
	* through get_control_int, which prints an error if it was never set.
	*
	* @return Control integer for imp_iz_recomb, or -1 if it was never set.
	*/
	int Options::imp_iz_recomb_int() const
		{return get_control_int("imp_iz_recomb", m_imp_iz_recomb_int);}

	/**
	* Returns the control integer for the tbound_type option.
	*
	* The control integer is an integer version of tbound_type that is
	* faster to check than the string: 0 = "absorbing", 1 = "periodic". It
	* is fetched through get_control_int, which prints an error if it was
	* never set.
	*
	* @return Control integer for tbound_type, or -1 if it was never set.
	*/
	int Options::tbound_type_int() const
		{return get_control_int("tbound_type", m_tbound_type_int);}

	/**
	* Returns the control integer for the min_xbound_type option.
	*
	* The control integer is an integer version of min_xbound_type that is
	* faster to check than the string: 0 = "absorbing", 1 = "periodic", 2 =
	* "core". It is fetched through get_control_int, which prints an error
	* if it was never set.
	*
	* @return Control integer for min_xbound_type, or -1 if it was never
	*         set.
	*/
	int Options::min_xbound_type_int() const
		{return get_control_int("min_xbound_type", m_min_xbound_type_int);}

	/**
	* Returns the control integer for the max_xbound_type option.
	*
	* The control integer is an integer version of max_xbound_type that is
	* faster to check than the string: 0 = "absorbing", 1 = "periodic", 2 =
	* "core". It is fetched through get_control_int, which prints an error
	* if it was never set.
	*
	* @return Control integer for max_xbound_type, or -1 if it was never
	*         set.
	*/
	int Options::max_xbound_type_int() const
		{return get_control_int("max_xbound_type", m_max_xbound_type_int);}

	/**
	* Returns the control integer for the min_ybound_type option.
	*
	* The control integer is an integer version of min_ybound_type that is
	* faster to check than the string: 0 = "absorbing", 1 = "periodic", 2 =
	* "core". It is fetched through get_control_int, which prints an error
	* if it was never set.
	*
	* @return Control integer for min_ybound_type, or -1 if it was never
	*         set.
	*/
	int Options::min_ybound_type_int() const
		{return get_control_int("min_ybound_type", m_min_ybound_type_int);}

	/**
	* Returns the control integer for the max_ybound_type option.
	*
	* The control integer is an integer version of max_ybound_type that is
	* faster to check than the string: 0 = "absorbing", 1 = "periodic", 2 =
	* "core". It is fetched through get_control_int, which prints an error
	* if it was never set.
	*
	* @return Control integer for max_ybound_type, or -1 if it was never
	*         set.
	*/
	int Options::max_ybound_type_int() const
		{return get_control_int("max_ybound_type", m_max_ybound_type_int);}

	/**
	* Returns the control integer for the min_zbound_type option.
	*
	* The control integer is an integer version of min_zbound_type that is
	* faster to check than the string: 0 = "absorbing", 1 = "periodic", 2 =
	* "core". It is fetched through get_control_int, which prints an error
	* if it was never set.
	*
	* @return Control integer for min_zbound_type, or -1 if it was never
	*         set.
	*/
	int Options::min_zbound_type_int() const
		{return get_control_int("min_zbound_type", m_min_zbound_type_int);}

	/**
	* Returns the control integer for the max_zbound_type option.
	*
	* The control integer is an integer version of max_zbound_type that is
	* faster to check than the string: 0 = "absorbing", 1 = "periodic", 2 =
	* "core". It is fetched through get_control_int, which prints an error
	* if it was never set.
	*
	* @return Control integer for max_zbound_type, or -1 if it was never
	*         set.
	*/
	int Options::max_zbound_type_int() const
		{return get_control_int("max_zbound_type", m_max_zbound_type_int);}

	/**
	* Returns the control integer for the calc_grad_elec option.
	*
	* The control integer is an integer version of calc_grad_elec that is
	* faster to check than the string: 0 = "off", 1 = "on". It is fetched
	* through get_control_int, which prints an error if it was never set.
	*
	* @return Control integer for calc_grad_elec, or -1 if it was never set.
	*/
	int Options::calc_grad_elec_int() const
		{return get_control_int("calc_grad_elec", m_calc_grad_elec_int);}

	// Setter declarations for internal control variables

	/**
	* Sets the control integer for the var_red_split option directly.
	*
	* Unlike the other control integers, this one can be changed without
	* changing the var_red_split string option.
	*
	* @param var_red_split_int New value of the control integer.
	*/
	void Options::set_var_red_split_int(int var_red_split_int)
		{m_var_red_split_int = var_red_split_int;}
	
#ifndef USE_CUDA

	// If these definitions gets called, it means you're in the CPU-only version
	// of the code trying to call GPU-specific code.
	// Defined in src/cuda/options.cu since it calls CUDA code

	/**
	* Placeholder for copying the options to a GPU in CPU-only builds.
	*
	* The real definition is in src/cuda/options.cu. This version only
	* prints an error saying GPU support was not compiled in.
	*
	* @param device_id ID of the GPU to copy to (unused).
	*/
	OptionsDevice* Options::to_device(int device_id)
	{
		std::cerr << "Error! Options::to_device was called but GPU support"
				  << " was not compiled in.\n";
	}

	/**
	* Placeholder for freeing device-side options in CPU-only builds.
	*
	* The real definition is in src/cuda/options.cu. This version only
	* prints an error saying GPU support was not compiled in.
	*
	* @param opts_d Device-side options to free (unused).
	* @param device_id ID of the GPU the options are on (unused).
	*/
	void free_opts(OptionsDevice* opts_d, int device_id)
	{
		std::cerr << "Error! Options::free_bkg was called but GPU support"
				  << " was not compiled in.\n";
	}

#endif

	/**
	* Checks a value against a list of valid values for an option.
	*
	* If the value is not in the list, an error naming the option and the
	* value is printed.
	*
	* @tparam T Type of the option value.
	* @param var Name of the option (for the error message).
	* @param value Value to check.
	* @param valid_values List of valid values, e.g., {"on", "off"}.
	* @return True if the value is valid, false otherwise.
	*/
	template <typename T>
	bool Options::check_input(const std::string& var, const T& value, 
		std::initializer_list<T> valid_values)
	{
		// Check input value against allowed values
		for (const auto& v : valid_values)
		{
			if (value == v) 
			{
				return true;
			}
		}

		// If allowed value not found, print error
		std::cerr << "Error! Value for " << var << " not recognized: "
			<< value << '\n';
		return false;
	}
}
