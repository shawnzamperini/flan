/**
* @file options.cpp
*/
#include <iostream>
#include <string>
#include <tuple>

#include "options.h"


namespace Options
{
	// Default mapc2p that assumes computational coordinates are already the
	// Cartesian coordinates.
	std::tuple<double, double, double> mapc2p_default 
		(double xc, double yc, double zc)
	{
		// Default implementation assumes the computational coordinates are
		// already Cartesian
		return std::make_tuple(xc, yc, zc);
	}

	// Constructors
	Options::Options()
	{}
	
	// Setter definitions. They're all the same, so just throwing together
	// as a mass of code. Can organize if necessary, but it would probably
	// just take up three times as much space.
	void Options::set_case_name(std::string case_name) 
		{m_case_name = case_name;}
	void Options::set_gkyl_dir(std::string gkyl_dir) 
		{m_gkyl_dir = gkyl_dir;}
	void Options::set_gkyl_casename(std::string gkyl_casename) 
		{m_gkyl_casename = gkyl_casename;}
	void Options::set_gkyl_frame_start(int gkyl_frame_start) 
		{m_gkyl_frame_start = gkyl_frame_start;}
	void Options::set_gkyl_frame_end(int gkyl_frame_end) 
		{m_gkyl_frame_end = gkyl_frame_end;}
	void Options::set_gkyl_elec_name(std::string gkyl_elec_name) 
		{m_gkyl_elec_name = gkyl_elec_name;}
	void Options::set_gkyl_ion_name(std::string gkyl_ion_name) 
		{m_gkyl_ion_name = gkyl_ion_name;}
	void Options::set_gkyl_elec_mass_amu(double gkyl_elec_mass_amu) 
		{m_gkyl_elec_mass_amu = gkyl_elec_mass_amu;}
	void Options::set_gkyl_ion_mass_amu(double gkyl_ion_mass_amu) 
		{m_gkyl_ion_mass_amu = gkyl_ion_mass_amu;}
	void Options::set_gkyl_file_type(std::string gkyl_file_type) 
		{m_gkyl_file_type = gkyl_file_type;}
	void Options::set_geotype(std::string geotype) 
		{
			m_geotype = geotype;
			if (geotype == "simple helical") m_geotype_int = 0;
			else if (geotype == "toroidal") m_geotype_int = 1;
			else
			{
				std::cerr << "Error! geotype not recognized: " << geotype 
					<< '\n';
			}
		}
	void Options::set_lcfs_x(double lcfs_x) 
		{m_lcfs_x = lcfs_x;}
	void Options::set_imp_atom_num(int imp_atom_num) 
		{m_imp_atom_num = imp_atom_num;}
	void Options::set_imp_mass_amu(double imp_mass_amu) 
		{m_imp_mass_amu = imp_mass_amu;}
	void Options::set_imp_init_charge(int imp_init_charge) 
		{m_imp_init_charge = imp_init_charge;}
	void Options::set_imp_num(int imp_num) 
		{m_imp_num = imp_num;}
	void Options::set_imp_xmin(double imp_xmin) 
		{m_imp_xmin = imp_xmin;}
	void Options::set_imp_xmax(double imp_xmax) 
		{m_imp_xmax = imp_xmax;}
	void Options::set_imp_ystart_opt(std::string imp_ystart_opt) 
		{
			m_imp_ystart_opt = imp_ystart_opt;
			if (imp_ystart_opt == "single_value") m_imp_ystart_opt_int = 0;
			else if (imp_ystart_opt == "range") m_imp_ystart_opt_int = 1;
			else
			{
				std::cerr << "Error! imp_ystart_opt not recognized: " 
					<< imp_ystart_opt << '\n';
			}
		}
	void Options::set_imp_ystart_val(double imp_ystart_val) 
		{m_imp_ystart_val = imp_ystart_val;}
	void Options::set_imp_zstart_opt(std::string imp_zstart_opt) 
		{
			m_imp_zstart_opt = imp_zstart_opt;
			if (imp_zstart_opt == "single_value") m_imp_zstart_opt_int = 0;
			else if (imp_zstart_opt == "range") m_imp_zstart_opt_int = 1;
		}
	void Options::set_imp_zstart_val(double imp_zstart_val) 
		{m_imp_zstart_val = imp_zstart_val;}
	void Options::set_imp_collisions(std::string imp_collisions) 
		{
			m_imp_collisions = imp_collisions;
			if (imp_collisions == "off") m_imp_collisions_int = 0;
			else if (imp_collisions == "on") m_imp_collisions_int = 1;
		}
	void Options::set_var_red(std::string var_red) 
		{
			m_var_red = var_red;
			if (var_red == "off") m_var_red_int = 0;
			else if (var_red == "on") m_var_red_int = 1;
		}

	void Options::set_var_red_freq(double var_red_freq) 
		{m_var_red_freq = var_red_freq;}
	
	void Options::set_var_red_min_weight(double var_red_min_weight) 
		{m_var_red_min_weight = var_red_min_weight;}
	
	void Options::set_var_red_med_mod(double var_red_med_mod) 
		{m_var_red_med_mod = var_red_med_mod;}
	
	void Options::set_imp_time_step_opt(std::string imp_time_step_opt) 
		{
			m_imp_time_step_opt = imp_time_step_opt;
			if (imp_time_step_opt == "constant") m_imp_time_step_opt_int = 0;
			else if (imp_time_step_opt == "variable") 
				m_imp_time_step_opt_int = 1;
		}
	
	void Options::set_imp_time_step(double imp_time_step) 
		{m_imp_time_step = imp_time_step;}
	void Options::set_imp_time_step_min(double imp_time_step_min) 
		{m_imp_time_step_min = imp_time_step_min;}
	void Options::set_imp_source_scale_fact(double imp_source_scale_fact) 
		{m_imp_source_scale_fact = imp_source_scale_fact;}
	void Options::set_imp_vel_stats(std::string imp_vel_stats) 
		{
			m_imp_vel_stats = imp_vel_stats;
			if (imp_vel_stats == "off") m_imp_vel_stats_int = 0;
			else if (imp_vel_stats == "on") m_imp_vel_stats_int = 1;
		}
	void Options::set_imp_xbound_buffer(double imp_xbound_buffer) 
		{m_imp_xbound_buffer = imp_xbound_buffer;}
	void Options::set_imp_iz_recomb(std::string imp_iz_recomb) 
		{
			m_imp_iz_recomb = imp_iz_recomb;
			if (imp_iz_recomb == "off") m_imp_iz_recomb_int = 0;
			else if (imp_iz_recomb == "on") m_imp_iz_recomb_int = 1;
		}
	void Options::set_print_interval(int print_interval) 
		{m_print_interval = print_interval;}
	void Options::set_openadas_root(std::string openadas_root) 
		{m_openadas_root = openadas_root;}
	void Options::set_openadas_year(int openadas_year) 
		{m_openadas_year = openadas_year;}
	void Options::set_mapc2p(Mapc2p_ptr mapc2p) 
		{m_mapc2p = mapc2p;}

	// Accessor definitions
	const std::string& Options::case_name() const 
		{return m_case_name;}
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
	const std::string& Options::geotype() const 
		{return m_geotype;}
	const double Options::lcfs_x() const 
		{return m_lcfs_x;}
	const int Options::imp_atom_num() const 
		{return m_imp_atom_num;}
	const double Options::imp_mass_amu() const 
		{return m_imp_mass_amu;}
	const int Options::imp_init_charge() const 
		{return m_imp_init_charge;}
	const int Options::imp_num() const
		{return m_imp_num;}
	const double Options::imp_xmin() const 
		{return m_imp_xmin;}
	const double Options::imp_xmax() const 
		{return m_imp_xmax;}
	const std::string& Options::imp_ystart_opt() const 
		{return m_imp_ystart_opt;}
	const double Options::imp_ystart_val() const 
		{return m_imp_ystart_val;}
	const std::string& Options::imp_zstart_opt() const 
		{return m_imp_zstart_opt;}
	const double Options::imp_zstart_val() const 
		{return m_imp_zstart_val;}
	const std::string& Options::imp_collisions() const 
		{return m_imp_collisions;}
	const std::string& Options::var_red() const 
		{return m_var_red;}
	const double Options::var_red_freq() const 
		{return m_var_red_freq;}
	const double Options::var_red_min_weight() const 
		{return m_var_red_min_weight;}
	const double Options::var_red_med_mod() const 
		{return m_var_red_med_mod;}
	const std::string& Options::imp_time_step_opt() const 
		{return m_imp_time_step_opt;}
	const double Options::imp_time_step() const 
		{return m_imp_time_step;}
	const double Options::imp_time_step_min() const 
		{return m_imp_time_step_min;}
	const double Options::imp_source_scale_fact() const 
		{return m_imp_source_scale_fact;}
	const std::string& Options::imp_vel_stats() const 
		{return m_imp_vel_stats;}
	const double Options::imp_xbound_buffer() const 
		{return m_imp_xbound_buffer;}
	const std::string& Options::imp_iz_recomb() const 
		{return m_imp_iz_recomb;}
	const int Options::print_interval() const 
		{return m_print_interval;}
	const std::string& Options::openadas_root() const 
		{return m_openadas_root;}
	const int Options::openadas_year() const 
		{return m_openadas_year;}
	const Mapc2p_ptr Options::mapc2p() const 
		{return m_mapc2p;}
	
	// Accessors for internal control variables
	const int Options::imp_ystart_opt_int() const
		{return m_imp_ystart_opt_int;}
	const int Options::imp_zstart_opt_int() const
		{return m_imp_zstart_opt_int;}
	const int Options::imp_collisions_int() const
		{return m_imp_collisions_int;}
	const int Options::var_red_int() const
		{return m_var_red_int;}
	const int Options::imp_time_step_opt_int() const
		{return m_imp_time_step_opt_int;}
	const int Options::imp_vel_stats_int() const
		{return m_imp_vel_stats_int;}
	const int Options::imp_iz_recomb_int() const
		{return m_imp_iz_recomb_int;}
	const int Options::geotype_int() const
		{return m_geotype_int;}
	

}
