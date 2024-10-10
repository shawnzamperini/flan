/**
* @file read_input.h
* @brief Header file for read_input.cpp
*/

#ifndef READ_INPUT_H
#define READ_INPUT_H

#include <array>
#include <string>
#include <fstream>

#include "input_classes.h"

namespace Input
{
	
	// Enumeration containing all the input options. There are used to index
	// the options array defined in read_input.cpp. These are here in the 
	// header file so other files can include this file to access them.
	enum Names
	{
		gkyl_dir,
		gkyl_casename,
		gkyl_frame_start,
		gkyl_frame_end,
		gkyl_elec_name,
		gkyl_ion_name,
		gkyl_elec_mass_amu,
		gkyl_ion_mass_amu,
		gkyl_file_type,
		imp_atom_num,
		imp_mass_amu,
		imp_init_charge,
		imp_num,
		imp_xmin,
		imp_xmax,
		imp_zstart_opt,
		imp_zstart_val,
		imp_collisions,
		imp_var_reduct,
		imp_var_reduct_freq,
		imp_var_reduct_min_weight,
		imp_time_step,
		imp_source_scale_fact,
		openadas_root,
		openadas_year,
		max_input_opts,  // Always leave this at the end of Names
	};

	extern std::array<OptionBase*, max_input_opts> options;

	// Provide a string and return the corresponding enumerator in Names
	int str_to_name(std::string_view input_str);

	std::string_view get_opt_str(Names name);
	int get_opt_int(Names name);
	double get_opt_dbl(Names name);

	// Takes in a line from an input file as a stringstream object and returns
	// it as a normal string for further manipulation. 
	std::string get_str_from_line(std::stringstream& input_line_stream);

	// Using the option name (input_name) and the value pulled from the input
	// file (input_val), save this value into the main options array. Returns
	// true if successful, false if not.
	bool save_input_opt(std::string_view input_name, std::string_view input_val);

	// Take a line from an input file and parse it, saving the input into 
	// the options array if the line is not blank or a comment. 
	void parse_input_line(std::string_view input_line);

	// This goes through all the input options and changes the values from 
	// their defaults using the value included in the input file. All the input 
	// options can be found in the options array.
	void load_input_opts(std::ifstream& input_stream);
}

#endif
