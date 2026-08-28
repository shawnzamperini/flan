/**
* @file options.cpp
*/
#include <iostream>
#include <string>
#include <tuple>
#include <initializer_list>

#include "options.h"
#include "options_device.h"


namespace Options
{

	// Constructors
	Options::Options()
	{
		// Set internal control integers based on (at this point, default)
		// options. They'll get changed later if the option changes.
		initialize_control_ints();
	}
	
	// Initialize the control integers, which are generally used to avoid
	// having to compare against strings in the code for speed. 
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
	void Options::set_case_name(std::string case_name) 
		{m_case_name = case_name;}

	// use_gpu
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

	// seed
	void Options::set_seed(int seed) 
		{m_seed = seed;}

	// slot_cap
	void Options::set_slot_cap(int slot_cap) 
		{m_slot_cap = slot_cap;}

	// bkg_source
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

	// test_opt
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

	// save_track
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
	
	// gkyl_dir
	void Options::set_gkyl_dir(std::string gkyl_dir) 
		{m_gkyl_dir = gkyl_dir;}
	
	// gkyl_casename
	void Options::set_gkyl_casename(std::string gkyl_casename) 
		{m_gkyl_casename = gkyl_casename;}
	
	// gkyl_frame_start
	void Options::set_gkyl_frame_start(int gkyl_frame_start) 
		{m_gkyl_frame_start = gkyl_frame_start;}
	
	// gkyl_frame_end
	void Options::set_gkyl_frame_end(int gkyl_frame_end) 
		{m_gkyl_frame_end = gkyl_frame_end;}

	// gkyl_elec_name
	void Options::set_gkyl_elec_name(std::string gkyl_elec_name) 
		{m_gkyl_elec_name = gkyl_elec_name;}
	
	// gkyl_ion_name
	void Options::set_gkyl_ion_name(std::string gkyl_ion_name) 
		{m_gkyl_ion_name = gkyl_ion_name;}
	
	// gkyl_elec_mass_amu
	void Options::set_gkyl_elec_mass_amu(double gkyl_elec_mass_amu) 
		{m_gkyl_elec_mass_amu = gkyl_elec_mass_amu;}
	
	// gkyl_ion_mass_amu
	void Options::set_gkyl_ion_mass_amu(double gkyl_ion_mass_amu) 
		{m_gkyl_ion_mass_amu = gkyl_ion_mass_amu;}

	// gkyl_file_type
	void Options::set_gkyl_file_type(std::string gkyl_file_type) 
	{
		if (check_input<std::string>("gkyl_file_type", gkyl_file_type,
			{"binary"}))
		{
			m_gkyl_file_type = gkyl_file_type;
		}
	}

	// gkyl_file_type
	void Options::set_gkyl_moment_type(std::string gkyl_moment_type) 
	{
		if (check_input<std::string>("gkyl_moment_type", gkyl_moment_type,
			{"bimaxwellian", "maxwellian", "m0m1m2"}))
		{
			m_gkyl_moment_type = gkyl_moment_type;
		}
	}

	// calc_grad_elec
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

	// tbound_type
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

	// min_xbound_type
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

	// max_xbound_type
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

	// min_ybound_type
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

	// max_ybound_type
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

	// min_zbound_type
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

	// max_zbound_type
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

	// lcfs_x
	void Options::set_lcfs_x(double lcfs_x) 
		{m_lcfs_x = lcfs_x;}

	// sep_x_bc_xp_z1
	void Options::set_sep_x_bc_xp_z1(double sep_x_bc_xp_z1) 
		{m_sep_x_bc_xp_z1 = sep_x_bc_xp_z1;}

	// sep_x_bc_xp_z2
	void Options::set_sep_x_bc_xp_z2(double sep_x_bc_xp_z2) 
		{m_sep_x_bc_xp_z2 = sep_x_bc_xp_z2;}

	// imp_atom_num
	void Options::set_imp_atom_num(int imp_atom_num) 
		{m_imp_atom_num = imp_atom_num;}

	// imp_mass_amu
	void Options::set_imp_mass_amu(double imp_mass_amu) 
		{m_imp_mass_amu = imp_mass_amu;}
	
	// imp_init_charge
	void Options::set_imp_init_charge(int imp_init_charge) 
		{m_imp_init_charge = imp_init_charge;}
	
	// imp_num
	void Options::set_imp_num(int imp_num) 
		{m_imp_num = imp_num;}
	
	// imp_tstart_opt
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

	// imp_tstart_val
	void Options::set_imp_tstart_val(double imp_tstart_val) 
		{m_imp_tstart_val = imp_tstart_val;}
		
	// imp_trange_min
	void Options::set_imp_trange_min(double imp_trange_min) 
		{m_imp_trange_min = imp_trange_min;}
	
	// imp_trange_max
	void Options::set_imp_trange_max(double imp_trange_max) 
		{m_imp_trange_max = imp_trange_max;}

	// imp_xstart_opt
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

	// imp_xstart_val
	void Options::set_imp_xstart_val(double imp_xstart_val) 
		{m_imp_xstart_val = imp_xstart_val;}
		
	// imp_xrange_min
	void Options::set_imp_xrange_min(double imp_xrange_min) 
		{m_imp_xrange_min = imp_xrange_min;}
	
	// imp_xrange_max
	void Options::set_imp_xrange_max(double imp_xrange_max) 
		{m_imp_xrange_max = imp_xrange_max;}
	
	// imp_ystart_opt
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
	
	// imp_ystart_val
	void Options::set_imp_ystart_val(double imp_ystart_val) 
		{m_imp_ystart_val = imp_ystart_val;}

	// imp_yrange_min
	void Options::set_imp_yrange_min(double imp_yrange_min) 
		{m_imp_yrange_min = imp_yrange_min;}
	
	// imp_yrange_max
	void Options::set_imp_yrange_max(double imp_yrange_max) 
		{m_imp_yrange_max = imp_yrange_max;}
	
	// imp_zstart_opt
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
	
	// imp_zstart_val
	void Options::set_imp_zstart_val(double imp_zstart_val) 
		{m_imp_zstart_val = imp_zstart_val;}

	// imp_zrange_min
	void Options::set_imp_zrange_min(double imp_zrange_min) 
		{m_imp_zrange_min = imp_zrange_min;}
	
	// imp_zrange_max
	void Options::set_imp_zrange_max(double imp_zrange_max) 
		{m_imp_zrange_max = imp_zrange_max;}
	
	// imp_temp_start_opt
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
	
	// imp_temp_start_val
	void Options::set_imp_temp_start_val(double imp_temp_start_val) 
		{m_imp_temp_start_val = imp_temp_start_val;}

	// imp_collisions
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

	// var_red_split
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

	// var_red_import
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

	// var_red_freq
	void Options::set_var_red_freq(double var_red_freq) 
		{m_var_red_freq = var_red_freq;}
	
	// var_red_min_weight
	void Options::set_var_red_min_weight(double var_red_min_weight) 
		{m_var_red_min_weight = var_red_min_weight;}
	
	// var_red_med_mod
	void Options::set_var_red_med_mod(double var_red_med_mod) 
		{m_var_red_med_mod = var_red_med_mod;}
	
	// var_red_rusrol
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

	// var_red_rusrol_prob
	void Options::set_var_red_rusrol_prob(double var_red_rusrol_prob) 
		{m_var_red_rusrol_prob = var_red_rusrol_prob;}

	// imp_time_step_opt
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
	
	// imp_time_step
	void Options::set_imp_time_step(double imp_time_step) 
		{m_imp_time_step = imp_time_step;}
	
	// imp_time_step_min
	void Options::set_imp_time_step_min(double imp_time_step_min) 
		{m_imp_time_step_min = imp_time_step_min;}
	
	// imp_source_scale_fact
	void Options::set_imp_source_scale_fact(double imp_source_scale_fact) 
		{m_imp_source_scale_fact = imp_source_scale_fact;}
	
	// imp_xbound_buffer
	void Options::set_imp_xbound_buffer(double imp_xbound_buffer) 
		{m_imp_xbound_buffer = imp_xbound_buffer;}

	// imp_ybound_buffer
	void Options::set_imp_ybound_buffer(double imp_ybound_buffer) 
		{m_imp_ybound_buffer = imp_ybound_buffer;}

	// imp_zbound_buffer
	void Options::set_imp_zbound_buffer(double imp_zbound_buffer) 
		{m_imp_zbound_buffer = imp_zbound_buffer;}
	
	// imp_iz_recomb
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
	
	// print_interval
	void Options::set_print_interval(int print_interval) 
		{m_print_interval = print_interval;}
	
	// openadas_root
	void Options::set_openadas_root(std::string openadas_root) 
		{m_openadas_root = openadas_root;}
	
	// openadas_year
	void Options::set_openadas_year(int openadas_year) 
		{m_openadas_year = openadas_year;}
	
	// Accessor definitions
	const std::string& Options::case_name() const 
		{return m_case_name;}
	const std::string& Options::use_gpu() const 
		{return m_use_gpu;}
	const int Options::seed() const 
		{return m_seed;}
	const int Options::slot_cap() const 
		{return m_slot_cap;}
	const std::string& Options::bkg_source() const 
		{return m_bkg_source;}
	const std::string& Options::test_opt() const 
		{return m_test_opt;}
	const std::string& Options::save_track() const 
		{return m_save_track;}
	const std::string& Options::gkyl_dir() const 
		{return m_gkyl_dir;}
	const std::string& Options::gkyl_casename() const 
		{return m_gkyl_casename;}
	const int Options::gkyl_frame_start() const 
		{return m_gkyl_frame_start;}
	const int Options::gkyl_frame_end() const 
		{return m_gkyl_frame_end;}
	const std::string& Options::gkyl_elec_name() const 
		{return m_gkyl_elec_name;}
	const std::string& Options::gkyl_ion_name() const 
		{return m_gkyl_ion_name;}
	const double Options::gkyl_elec_mass_amu() const 
		{return m_gkyl_elec_mass_amu;}
	const double Options::gkyl_ion_mass_amu() const 
		{return m_gkyl_ion_mass_amu;}
	const std::string& Options::gkyl_file_type() const 
		{return m_gkyl_file_type;}
	const std::string& Options::gkyl_moment_type() const 
		{return m_gkyl_moment_type;}
	const std::string& Options::min_xbound_type() const 
		{return m_min_xbound_type;}
	const std::string& Options::max_xbound_type() const 
		{return m_max_xbound_type;}
	const std::string& Options::min_ybound_type() const 
		{return m_min_ybound_type;}
	const std::string& Options::max_ybound_type() const 
		{return m_max_ybound_type;}
	const std::string& Options::min_zbound_type() const 
		{return m_min_zbound_type;}
	const std::string& Options::max_zbound_type() const 
		{return m_max_zbound_type;}
	const double Options::lcfs_x() const 
		{return m_lcfs_x;}
	const double Options::sep_x_bc_xp_z1() const 
		{return m_sep_x_bc_xp_z1;}
	const double Options::sep_x_bc_xp_z2() const 
		{return m_sep_x_bc_xp_z2;}
	const double Options::imp_xbound_buffer() const 
		{return m_imp_xbound_buffer;}
	const double Options::imp_ybound_buffer() const 
		{return m_imp_ybound_buffer;}
	const double Options::imp_zbound_buffer() const 
		{return m_imp_zbound_buffer;}
	const int Options::imp_atom_num() const 
		{return m_imp_atom_num;}
	const double Options::imp_mass_amu() const 
		{return m_imp_mass_amu;}
	const int Options::imp_init_charge() const 
		{return m_imp_init_charge;}
	const int Options::imp_num() const
		{return m_imp_num;}
	const std::string& Options::imp_tstart_opt() const 
		{return m_imp_tstart_opt;}
	const double Options::imp_tstart_val() const 
		{return m_imp_tstart_val;}
	const double Options::imp_trange_min() const 
		{return m_imp_trange_min;}
	const double Options::imp_trange_max() const 
		{return m_imp_trange_max;}
	const std::string& Options::imp_xstart_opt() const 
		{return m_imp_xstart_opt;}
	const double Options::imp_xstart_val() const 
		{return m_imp_xstart_val;}
	const double Options::imp_xrange_min() const 
		{return m_imp_xrange_min;}
	const double Options::imp_xrange_max() const 
		{return m_imp_xrange_max;}
	const std::string& Options::imp_ystart_opt() const 
		{return m_imp_ystart_opt;}
	const double Options::imp_ystart_val() const 
		{return m_imp_ystart_val;}
	const double Options::imp_yrange_min() const 
		{return m_imp_yrange_min;}
	const double Options::imp_yrange_max() const 
		{return m_imp_yrange_max;}
	const std::string& Options::imp_zstart_opt() const 
		{return m_imp_zstart_opt;}
	const double Options::imp_zstart_val() const 
		{return m_imp_zstart_val;}
	const double Options::imp_zrange_min() const 
		{return m_imp_zrange_min;}
	const double Options::imp_zrange_max() const 
		{return m_imp_zrange_max;}
	const std::string& Options::imp_temp_start_opt() const 
		{return m_imp_temp_start_opt;}
	const double Options::imp_temp_start_val() const 
		{return m_imp_temp_start_val;}
	const std::string& Options::imp_collisions() const 
		{return m_imp_collisions;}
	const std::string& Options::var_red_split() const 
		{return m_var_red_split;}
	const std::string& Options::var_red_import() const 
		{return m_var_red_import;}
	const double Options::var_red_freq() const 
		{return m_var_red_freq;}
	const double Options::var_red_min_weight() const 
		{return m_var_red_min_weight;}
	const double Options::var_red_med_mod() const 
		{return m_var_red_med_mod;}
	const std::string& Options::var_red_rusrol() const 
		{return m_var_red_rusrol;}
	const double Options::var_red_rusrol_prob() const 
		{return m_var_red_rusrol_prob;}
	const std::string& Options::imp_time_step_opt() const 
		{return m_imp_time_step_opt;}
	const double Options::imp_time_step() const 
		{return m_imp_time_step;}
	const double Options::imp_time_step_min() const 
		{return m_imp_time_step_min;}
	const double Options::imp_source_scale_fact() const 
		{return m_imp_source_scale_fact;}
	const std::string& Options::imp_iz_recomb() const 
		{return m_imp_iz_recomb;}
	const int Options::print_interval() const 
		{return m_print_interval;}
	const std::string& Options::openadas_root() const 
		{return m_openadas_root;}
	const int Options::openadas_year() const 
		{return m_openadas_year;}
	
	// Wrapper for returning control integer value that ensures they were
	// set to the respective option (i.e., not -1).
	const int Options::get_control_int(const std::string& var_name, 
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
	const int Options::use_gpu_int() const
		{return get_control_int("use_gpu", m_use_gpu_int);}
	const int Options::bkg_source_int() const
		{return get_control_int("bkg_source", m_bkg_source_int);}
	const int Options::test_opt_int() const
		{return get_control_int("test_opt", m_test_opt_int);}
	const int Options::save_track_int() const
		{return get_control_int("save_track", m_save_track_int);}
	const int Options::imp_tstart_opt_int() const
		{return get_control_int("imp_tstart_opt", m_imp_tstart_opt_int);}
	const int Options::imp_xstart_opt_int() const
		{return get_control_int("imp_xstart_opt", m_imp_xstart_opt_int);}
	const int Options::imp_ystart_opt_int() const
		{return get_control_int("imp_ystart_opt", m_imp_ystart_opt_int);}
	const int Options::imp_zstart_opt_int() const
		{return get_control_int("imp_zstart_opt", m_imp_zstart_opt_int);}
	const int Options::imp_temp_start_opt_int() const
		{return get_control_int("imp_temp_start_opt", m_imp_temp_start_opt_int);}
	const int Options::imp_collisions_int() const
		{return get_control_int("imp_collisions", m_imp_collisions_int);}
	const int Options::var_red_split_int() const
		{return get_control_int("var_red_split", m_var_red_split_int);}
	const int Options::var_red_import_int() const
		{return get_control_int("var_red_import", m_var_red_import_int);}
	const int Options::var_red_rusrol_int() const
		{return get_control_int("var_red_rusrol", m_var_red_rusrol_int);}
	const int Options::imp_time_step_opt_int() const
		{return get_control_int("imp_time_step_opt", m_imp_time_step_opt_int);}
	const int Options::imp_iz_recomb_int() const
		{return get_control_int("imp_iz_recomb", m_imp_iz_recomb_int);}
	const int Options::tbound_type_int() const
		{return get_control_int("tbound_type", m_tbound_type_int);}
	const int Options::min_xbound_type_int() const
		{return get_control_int("min_xbound_type", m_min_xbound_type_int);}
	const int Options::max_xbound_type_int() const
		{return get_control_int("max_xbound_type", m_max_xbound_type_int);}
	const int Options::min_ybound_type_int() const
		{return get_control_int("min_ybound_type", m_min_ybound_type_int);}
	const int Options::max_ybound_type_int() const
		{return get_control_int("max_ybound_type", m_max_ybound_type_int);}
	const int Options::min_zbound_type_int() const
		{return get_control_int("min_zbound_type", m_min_zbound_type_int);}
	const int Options::max_zbound_type_int() const
		{return get_control_int("max_zbound_type", m_max_zbound_type_int);}
	const int Options::calc_grad_elec_int() const
		{return get_control_int("calc_grad_elec", m_calc_grad_elec_int);}

	// Setter declarations for internal control variables
	void Options::set_var_red_split_int(int var_red_split_int)
		{m_var_red_split_int = var_red_split_int;}
	
#ifndef USE_CUDA

	// If these definitions gets called, it means you're in the CPU-only version
	// of the code trying to call GPU-specific code.
	// Defined in src/cuda/options.cu since it calls CUDA code
	OptionsDevice* Options::to_device(int device_id)
	{
		std::cerr << "Error! Options::to_device was called but GPU support"
				  << " was not compiled in.\n";
	}

	void free_opts(OptionsDevice* opts_d, int device_id)
	{
		std::cerr << "Error! Options::free_bkg was called but GPU support"
				  << " was not compiled in.\n";
	}

#endif

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
