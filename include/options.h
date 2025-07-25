#ifndef OPTIONS_H
#define OPTIONS_H

/**
* @file options.h
*/
#include <tuple>

#include "flan_types.h"


namespace Options
{
	/** 
	* @brief Default mapc2p that assumes computational coordinates are 
	*        already the Cartesian coordinates.
	*/
	std::tuple<double, double, double> mapc2p_default 
		(double xc, double yc, double zc);

	/**
	* @brief Class for storing all the options that control a simulation.
	*
	* Default values for every option are included in the definitions of each
	* member variable. These can then be overwritten as in @ref 
	* load_options.cpp
	*/
	class Options
	{
	private:

		// Input options and their defaults are defined here

		// Meta options for the simulation
		std::string m_case_name              {"undefined"};

		// Options related to reading in Gkeyll files
		std::string m_gkyl_dir               {"undefined"};
		std::string m_gkyl_casename          {"undefined"};
		int m_gkyl_frame_start                         {0};
		int m_gkyl_frame_end                           {1};
		std::string m_gkyl_elec_name               {"elc"};
		std::string m_gkyl_ion_name                {"ion"};
		double m_gkyl_elec_mass_amu             {0.000548};
		double m_gkyl_ion_mass_amu                 {2.014};
		std::string m_gkyl_file_type            {"binary"};
		std::string m_gkyl_moment_type    {"bimaxwellian"};

		// Geometry options
		double m_lcfs_x                              {0.0};
		double m_imp_xbound_buffer                   {0.0};
		std::string m_min_xbound_type        {"absorbing"};

		// Impurity chracteristics
		int m_imp_atom_num                            {74};
		double m_imp_mass_amu                     {183.84};
		int m_imp_init_charge                          {1};

		// Impurity transport options
		int m_imp_num                                  {1};
		double m_imp_xmin                            {0.0};
		double m_imp_xmax                            {0.0};
		std::string m_imp_ystart_opt      {"single_value"};
		double m_imp_ystart_val                      {0.0};
		std::string m_imp_zstart_opt      {"single_value"};
		double m_imp_zstart_val                      {0.0};
		std::string m_imp_collisions               {"off"};
		std::string m_imp_time_step_opt       {"variable"};
		double m_imp_time_step                      {1e-7};
		double m_imp_time_step_min                 {1e-12};
		double m_imp_source_scale_fact               {1.0};
		std::string m_imp_vel_stats                {"off"};
		std::string m_imp_iz_recomb                 {"on"};
		int m_print_interval                          {10};

		// Variance reduction options
		std::string m_var_red_split                {"off"};
		std::string m_var_red_import            {"median"};
		double m_var_red_freq                        {0.1};
		double m_var_red_min_weight                  {0.1};
		double m_var_red_med_mod                     {1.0};
		std::string m_var_red_rusrol               {"off"};
		double m_var_red_rusrol_prob				 {0.5};

		// OpenADAS options
		std::string m_openadas_root          {"undefined"};
		int m_openadas_year                           {50};

		// Function to map between computational and physical (Cartesian)
		// coordinates. Generally copy/pasted from Gkeyll input file. 
		// mapc2p_default is defined at the beginning of this file before the
		// class. Mapc2p_ptr is defined in flan_types.h.
		Mapc2p_ptr m_mapc2p {&mapc2p_default};

		// Internal control variables for string options. The string options
		// are just there to be user-friendly, but it is more efficient to
		// use integers instead of comparing strings all the time. These are
		// all set automatically within the corresponding setters.
		int m_imp_ystart_opt_int	{0};
		int m_imp_zstart_opt_int	{0};
		int m_imp_collisions_int	{0};
		int m_var_red_split_int		{0};
		int m_var_red_import_int	{0};
		int m_var_red_rusrol_int	{0};
		int m_imp_time_step_opt_int {0};
		int m_imp_vel_stats_int		{0};
		int m_imp_iz_recomb_int		{0};
		int m_geotype_int			{0};
		int m_min_xbound_type_int	{0};

		/**
		* @brief Check input against a list of valid options
		*
		* @param var String containing variable name
		* @param value Template parameter of variable value
		* @param valid_values initializer_list of valid values, 
		*        e.g., {"on", "off"}
		*
		* @return true if valid, false if not
		*/
		template <typename T>
		bool check_input(const std::string& var, const T& value, 
			std::initializer_list<T> valid_values);

	public:

		Options();

		// Setter declarations
		void set_case_name(std::string case_name);
		void set_gkyl_dir(std::string gkyl_dir);
		void set_gkyl_casename(std::string gkyl_casename);
		void set_gkyl_frame_start(int gkyl_frame_start);
		void set_gkyl_frame_end(int gkyl_frame_end);
		void set_gkyl_elec_name(std::string gkyl_elec_name);
		void set_gkyl_ion_name(std::string gkyl_ion_name);
		void set_gkyl_elec_mass_amu(double gkyl_elec_mass_amu);
		void set_gkyl_ion_mass_amu(double gkyl_ion_mass_amu);
		void set_gkyl_file_type(std::string gkyl_file_type);
		void set_gkyl_moment_type(std::string gkyl_moment_type);
		void set_min_xbound_type(std::string min_xbound_type);
		void set_lcfs_x(double lcfs_x);
		void set_imp_atom_num(int imp_atom_num);
		void set_imp_mass_amu(double imp_mass_amu);
		void set_imp_init_charge(int imp_init_charge);
		void set_imp_num(int imp_num);
		void set_imp_xmin(double imp_xmin);
		void set_imp_xmax(double imp_xmax);
		void set_imp_ystart_opt(std::string imp_ystart_opt);
		void set_imp_ystart_val(double imp_ystart_val);
		void set_imp_zstart_opt(std::string imp_zstart_opt);
		void set_imp_zstart_val(double imp_zstart_val);
		void set_imp_collisions(std::string imp_collisions);
		void set_var_red_split(std::string var_red_split);
		void set_var_red_import(std::string var_red_import);
		void set_var_red_freq(double var_red_freq);
		void set_var_red_min_weight(double var_red_min_weight);
		void set_var_red_med_mod(double var_red_med_mod);
		void set_var_red_rusrol(std::string var_red_rusrol);
		void set_var_red_rusrol_prob(double var_red_rusrol_prob);
		void set_imp_time_step_opt(std::string imp_time_step_opt);
		void set_imp_time_step(double imp_time_step);
		void set_imp_time_step_min(double imp_time_step_min);
		void set_imp_source_scale_fact(double imp_source_scale_fact);
		void set_imp_vel_stats(std::string imp_vel_stats);
		void set_imp_xbound_buffer(double imp_xbound_buffer);
		void set_imp_iz_recomb(std::string imp_iz_recomb);
		void set_print_interval(int print_interval);
		void set_openadas_root(std::string openadas_root);
		void set_openadas_year(int openadas_year);
		void set_mapc2p(Mapc2p_ptr mapc2p);

		// Accessor declarations
		const std::string& case_name() const;
		const std::string& gkyl_dir() const;
		const std::string& gkyl_casename() const;
		const int gkyl_frame_start() const;
		const int gkyl_frame_end() const;
		const std::string& gkyl_elec_name() const;
		const std::string& gkyl_ion_name() const;
		const double gkyl_elec_mass_amu() const;
		const double gkyl_ion_mass_amu() const;
		const std::string& gkyl_file_type() const;
		const std::string& gkyl_moment_type() const;
		const std::string& min_xbound_type() const;
		const double lcfs_x() const;
		const int imp_atom_num() const;
		const double imp_mass_amu() const;
		const int imp_init_charge() const;
		const int imp_num() const;
		const double imp_xmin() const;
		const double imp_xmax() const;
		const std::string& imp_ystart_opt() const;
		const double imp_ystart_val() const;
		const std::string& imp_zstart_opt() const;
		const double imp_zstart_val() const;
		const std::string& imp_collisions() const;
		const std::string& var_red_split() const;
		const std::string& var_red_import() const;
		const double var_red_freq() const;
		const double var_red_min_weight() const;
		const double var_red_med_mod() const;
		const std::string& var_red_rusrol() const;
		const double var_red_rusrol_prob() const;
		const std::string& imp_time_step_opt() const;
		const double imp_time_step() const;
		const double imp_time_step_min() const;
		const double imp_source_scale_fact() const;
		const std::string& imp_vel_stats() const;
		const double imp_xbound_buffer() const;
		const std::string& imp_iz_recomb() const;
		const int print_interval() const;
		const std::string& openadas_root() const;
		const int openadas_year() const;
		const Mapc2p_ptr mapc2p() const;

		// Accessor declarations for internal control variables
		const int imp_ystart_opt_int() const;
		const int imp_zstart_opt_int() const;
		const int imp_collisions_int() const;
		const int var_red_split_int() const;
		const int var_red_import_int() const;
		const int var_red_rusrol_int() const;
		const int imp_time_step_opt_int() const;
		const int imp_vel_stats_int() const;
		const int imp_iz_recomb_int() const;
		const int geotype_int() const;
		const int min_xbound_type_int() const;

		// Setter declarations for internal control variables. Added as needed.
		void set_var_red_split_int(int var_red_split_int);
	};

}

#endif
