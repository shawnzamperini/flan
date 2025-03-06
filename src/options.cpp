/**
* @file options.cpp
*/
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
		{m_imp_ystart_opt = imp_ystart_opt;}
	void Options::set_imp_ystart_val(double imp_ystart_val) 
		{m_imp_ystart_val = imp_ystart_val;}
	void Options::set_imp_zstart_opt(std::string imp_zstart_opt) 
		{m_imp_zstart_opt = imp_zstart_opt;}
	void Options::set_imp_zstart_val(double imp_zstart_val) 
		{m_imp_zstart_val = imp_zstart_val;}
	void Options::set_imp_collisions(std::string imp_collisions) 
		{m_imp_collisions = imp_collisions;}
	void Options::set_imp_var_reduct(std::string imp_var_reduct) 
		{m_imp_var_reduct = imp_var_reduct;}
	void Options::set_imp_var_reduct_freq(double imp_var_reduct_freq) 
		{m_imp_var_reduct_freq = imp_var_reduct_freq;}
	void Options::set_imp_var_reduct_min_weight(double imp_var_reduct_min_weight) 
		{m_imp_var_reduct_min_weight = imp_var_reduct_min_weight;}
	void Options::set_imp_time_step_opt(std::string imp_time_step_opt) 
		{m_imp_time_step_opt = imp_time_step_opt;}
	void Options::set_imp_time_step(double imp_time_step) 
		{m_imp_time_step = imp_time_step;}
	void Options::set_imp_time_step_min(double imp_time_step_min) 
		{m_imp_time_step_min = imp_time_step_min;}
	void Options::set_imp_source_scale_fact(double imp_source_scale_fact) 
		{m_imp_source_scale_fact = imp_source_scale_fact;}
	void Options::set_imp_vel_stats(std::string imp_vel_stats) 
		{m_imp_vel_stats = imp_vel_stats;}
	void Options::set_imp_xbound_buffer(double imp_xbound_buffer) 
		{m_imp_xbound_buffer = imp_xbound_buffer;}
	void Options::set_imp_iz_recomb(std::string imp_iz_recomb) 
		{m_imp_iz_recomb = imp_iz_recomb;}
	void Options::set_openadas_root(std::string openadas_root) 
		{m_openadas_root = openadas_root;}
	void Options::set_openadas_year(int openadas_year) 
		{m_openadas_year = openadas_year;}
	void Options::set_mapc2p(Mapc2p_ptr mapc2p) 
		{m_mapc2p = mapc2p;}

	// Accessor definitions
	std::string Options::gkyl_dir() 
		{return m_gkyl_dir;}
	int Options::gkyl_frame_start() 
		{return m_gkyl_frame_start;}
	int Options::gkyl_frame_end() 
		{return m_gkyl_frame_end;}
	std::string Options::gkyl_elec_name() 
		{return m_gkyl_elec_name;}
	std::string Options::gkyl_ion_name() 
		{return m_gkyl_ion_name;}
	double Options::gkyl_elec_mass_amu() 
		{return m_gkyl_elec_mass_amu;}
	double Options::gkyl_ion_mass_amu() 
		{return m_gkyl_ion_mass_amu;}
	std::string Options::gkyl_file_type() 
		{return m_gkyl_file_type;}
	int Options::imp_atom_num() 
		{return m_imp_atom_num;}
	double Options::imp_mass_amu() 
		{return m_imp_mass_amu;}
	int Options::imp_init_charge() 
		{return m_imp_init_charge;}
	int Options::imp_num()
		{return m_imp_num;}
	double Options::imp_xmin() 
		{return m_imp_xmin;}
	double Options::imp_xmax() 
		{return m_imp_xmax;}
	std::string Options::imp_ystart_opt() 
		{return m_imp_ystart_opt;}
	double Options::imp_ystart_val() 
		{return m_imp_ystart_val;}
	std::string Options::imp_zstart_opt() 
		{return m_imp_zstart_opt;}
	double Options::imp_zstart_val() 
		{return m_imp_zstart_val;}
	std::string Options::imp_collisions() 
		{return m_imp_collisions;}
	std::string Options::imp_var_reduct() 
		{return m_imp_var_reduct;}
	double Options::imp_var_reduct_freq() 
		{return m_imp_var_reduct_freq;}
	double Options::imp_var_reduct_min_weight() 
		{return m_imp_var_reduct_min_weight;}
	std::string Options::imp_time_step_opt() 
		{return m_imp_time_step_opt;}
	double Options::imp_time_step() 
		{return m_imp_time_step;}
	double Options::imp_time_step_min() 
		{return m_imp_time_step_min;}
	double Options::imp_source_scale_fact() 
		{return m_imp_source_scale_fact;}
	std::string Options::imp_vel_stats() 
		{return m_imp_vel_stats;}
	double Options::imp_xbound_buffer() 
		{return m_imp_xbound_buffer;}
	std::string Options::imp_iz_recomb() 
		{return m_imp_iz_recomb;}
	std::string Options::openadas_root() 
		{return m_openadas_root;}
	int Options::openadas_year() 
		{return m_openadas_year;}
	Mapc2p_ptr Options::mapc2p() 
		{return m_mapc2p;}
	

}
