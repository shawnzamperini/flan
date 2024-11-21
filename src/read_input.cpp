/**
* @file read_input.cpp
* @brief Definitions of the input options and the routine to read them 
*
* This contains the declarations of everything needed to read in input from
* an input file into Option classes used by Flan to keep track of things.
* An input file takes the form of, e.g.,
* 
* @code
* # This is a comment
* string_option_name  |  name     Note: Do not use quotation marks.
* int_option_name     |  3
* double_option_name  |  183.34
* @endcode
* 
* The options can go in any order, and you can move the | around to make 
* everything align if you want to. The options must be valid options though or
* it will yell at you! For a list of all the options and their defaults, 
* see the options array in below (or the documentation). 
*/

#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "input_classes.h"
#include "read_input.h"


namespace Input
{
	/**
	* @brief Array containing all the input options along with default values.
	*
	* Add new input options here. Make sure to add it as a variable in 
	* Input::Names located in read_input.cpp
	*/
	std::array<OptionBase*, max_input_opts> options 
	{
		new OptionStr{"gkyl_dir",            "undefined"},
		new OptionStr{"gkyl_casename",       "undefined"},
		new OptionInt{"gkyl_frame_start",              0},
		new OptionInt{"gkyl_frame_end",                1},
		new OptionStr{"gkyl_elec_name",            "elc"},
		new OptionStr{"gkyl_ion_name",             "ion"},
		new OptionDbl{"gkyl_elec_mass_amu",     0.000548},
		new OptionDbl{"gkyl_ion_mass_amu",         2.014},
		new OptionStr{"gkyl_file_type",         "binary"},
		new OptionInt{"imp_atom_num",                 74},
		new OptionDbl{"imp_mass_amu",             183.84},
		new OptionInt{"imp_init_charge",               1},
		new OptionInt{"imp_num",                       1},
		new OptionDbl{"imp_xmin",                    0.0},
		new OptionDbl{"imp_xmax",                    0.0},
		new OptionStr{"imp_ystart_opt",   "single_value"},
		new OptionDbl{"imp_ystart_val",              0.0},
		new OptionStr{"imp_zstart_opt",   "single_value"},
		new OptionDbl{"imp_zstart_val",              0.0},
		new OptionStr{"imp_collisions",            "off"},
		new OptionStr{"imp_var_reduct",            "off"},
		new OptionDbl{"imp_var_reduct_freq",         0.1},
		new OptionDbl{"imp_var_reduct_min_weight",   0.1},
		new OptionStr{"imp_time_step_opt",    "constant"},
		new OptionDbl{"imp_time_step",              1e-7},
		new OptionDbl{"imp_time_step_min",         1e-12},
		new OptionDbl{"imp_source_scale_fact",		 1.0},
		new OptionStr{"imp_vel_stats",             "off"},
		new OptionDbl{"imp_xbound_buffer",		     0.0},
		new OptionStr{"imp_iz_recomb",             "yes"},
		new OptionStr{"openadas_root",       "undefined"},
		new OptionInt{"openadas_year",                50}
	};

	/**
	* @brief Provide a string and return the corresponding enumerator in Names
	*
	* @param input_str String containing name of input variable
	* @return Return the corresponding int value that identifies the variable
	* in Input::Names (see read_input.h). This does NOT return the actual value
	* of the input option.
	*/
	int str_to_name(std::string_view input_str)
	{
		if (input_str == "gkyl_dir")                   return gkyl_dir;
		if (input_str == "gkyl_casename")              return gkyl_casename;
		if (input_str == "gkyl_frame_start")           return gkyl_frame_start;
		if (input_str == "gkyl_frame_end")             return gkyl_frame_end;
		if (input_str == "gkyl_elec_name")             return gkyl_elec_name;
		if (input_str == "gkyl_ion_name")              return gkyl_ion_name;
		if (input_str == "gkyl_elec_mass_amu")         return gkyl_elec_mass_amu;
		if (input_str == "gkyl_ion_mass_amu")          return gkyl_ion_mass_amu;
		if (input_str == "gkyl_file_type")             return gkyl_file_type;
		if (input_str == "imp_atom_num")               return imp_atom_num;
		if (input_str == "imp_mass_amu")               return imp_mass_amu;
		if (input_str == "imp_init_charge")            return imp_init_charge;
		if (input_str == "imp_num")                    return imp_num;
		if (input_str == "imp_xmin")                   return imp_xmin;
		if (input_str == "imp_xmax")                   return imp_xmax;
		if (input_str == "imp_ystart_opt")             return imp_ystart_opt;
		if (input_str == "imp_ystart_val")             return imp_ystart_val;
		if (input_str == "imp_zstart_opt")             return imp_zstart_opt;
		if (input_str == "imp_zstart_val")             return imp_zstart_val;
		if (input_str == "imp_collisions")             return imp_collisions;
		if (input_str == "imp_var_reduct")             return imp_var_reduct;
		if (input_str == "imp_var_reduct_freq")        return imp_var_reduct_freq;
		if (input_str == "imp_var_reduct_min_weight")  return imp_var_reduct_min_weight;
		if (input_str == "imp_time_step_opt")          return imp_time_step_opt;
		if (input_str == "imp_time_step")              return imp_time_step;
		if (input_str == "imp_time_step_min")          return imp_time_step_min;
		if (input_str == "imp_source_scale_fact")      return imp_source_scale_fact;
		if (input_str == "imp_vel_stats")              return imp_vel_stats;
		if (input_str == "imp_xbound_buffer")          return imp_xbound_buffer;
		if (input_str == "imp_iz_recomb")              return imp_iz_recomb;
		if (input_str == "openadas_root")              return openadas_root;
		if (input_str == "openadas_year")              return openadas_year;
		else return -1;
	}

	// Should be as many Names as there are options in the array.
	// Note: It is on the developer that the order of the array and the names
	// in Names are the same! The actual input file can be in whatever order. 
	static_assert(std::ssize(options) == max_input_opts);

	/**
	* @brief Get string option
	*
	* Helper function to provide a Name from the enumeration above and it will
	* return the value held by the option. One function for strings, integers
	* and doubles. 
	* The reason we do not use a template here is because the type to be 
	* returned can not be deduced with CTAD.
	*
	* @param name The value in Input::Names corresponding to the string input
	* option.
	* @return Returns the input option as a string_view. 
	*/
	std::string_view get_opt_str(Names name)
	{
		OptionBase* opt_base_ptr {options[name]};
		return dynamic_cast<OptionStr*>(opt_base_ptr)->get_value();
	}

	/**
	* @brief Get int option
	*
	* Helper function to provide a Name from the enumeration above and it will
	* return the value held by the option. One function for strings, integers
	* and doubles. 
	* The reason we do not use a template here is because the type to be 
	* returned can not be deduced with CTAD.
	*
	* @param name The value in Input::Names corresponding to the int input
	* option.
	* @return Returns the input option as a int. 
	*/
	int get_opt_int(Names name)
	{
		OptionBase* opt_base_ptr {options[name]};
		return dynamic_cast<OptionInt*>(opt_base_ptr)->get_value();
	}

	/**
	* @brief Get double option
	*
	* Helper function to provide a Name from the enumeration above and it will
	* return the value held by the option. One function for strings, integers
	* and doubles. 
	* The reason we do not use a template here is because the type to be 
	* returned can not be deduced with CTAD.
	*
	* @param name The value in Input::Names corresponding to the double input
	* option.
	* @return Returns the input option as a double. 
	*/
	double get_opt_dbl(Names name)
	{
		OptionBase* opt_base_ptr {options[name]};
		return dynamic_cast<OptionDbl*>(opt_base_ptr)->get_value();
	}

	/**
	* @brief Get entire line from input file as a string
	*
	* Takes in a line from an input file as a stringstream object and returns
	* it as a normal string for further manipulation. 
	*
	* @param input_line_stream The line as a stringstream object
	* @return Returns the line converted to a string object
	*/
	std::string get_str_from_line(std::stringstream& input_line_stream)
	{
		std::string str {};
		std::getline(input_line_stream, str, '|');
		std::string::iterator end_pos = std::remove(str.begin(), 
			str.end(), ' ');
		str.erase(end_pos, str.end());
		return str;
	}

	/**
	* @brief Save input option into options array, overwriting the default
	*
	* Using the option name (input_name) and the value pulled from the input
	* file (input_val), save this value into the main options array. Returns
	* true if successful, false if not.
	*
	* @param input_name The name of the input option as a string 
	* @param input_val The value of the input option as a string
	* @return Returns true if opration successful, false if not
	*/
	bool save_input_opt(std::string_view input_name, 
		std::string_view input_val)
	{
			// Assign each variable to its appropriate spot in options.
			//std::cout << Input::str_to_name(input_name) << '\n';
			Input::OptionBase* opt_base_ptr {Input::options[
				Input::str_to_name(input_name)]};

			// Attempt to downcast to each inherited class until we hit upon 
			// the right one. Once we do, call put_value to overwrite the 
			// default value.
			if (dynamic_cast<Input::OptionStr*>(opt_base_ptr))
			{
				//std::cout << "String option\n";
				dynamic_cast<Input::OptionStr*>(opt_base_ptr)->
					put_value(input_val);
				//std::cout << *dynamic_cast<Input::OptionStr*>(opt_base_ptr) 
				//	<< '\n';
			}
			else if (dynamic_cast<Input::OptionInt*>(opt_base_ptr))
			{
				//std::cout << "Int option\n";
				dynamic_cast<Input::OptionInt*>(opt_base_ptr)->
					put_value(input_val);
				//std::cout << *dynamic_cast<Input::OptionInt*>(opt_base_ptr) 
				//	<< '\n';
			}
			else if (dynamic_cast<Input::OptionDbl*>(opt_base_ptr))
			{
				//std::cout << "Double option\n";
				dynamic_cast<Input::OptionDbl*>(opt_base_ptr)->
					put_value(input_val);
				//std::cout << *dynamic_cast<Input::OptionDbl*>(opt_base_ptr) 
				//	<< '\n';
			}
			else
			{
				std::cout << "Error: could not save input value\n";
				return false;
			}

			return true;
	}

	/**
	* @brief Helper function to determine if string is a comment or blank
	*
	* Comments start with a #
	*
	* @return Return true if comment or blank, false otherwise
	*/
	bool not_comment_or_blank_line(const std::string_view line) 
	{
		return !line.starts_with("#") && 
			(line.find_first_not_of(' ') != std::string::npos);
	}

	/**
	* @brief Parse line from input file, saving input option if valid
	*
	* Take a line from an input file and parse it, saving the input into 
	* the options array if the line is not blank or a comment. 
	*
	* @param input_line Line from input file as a string_view
	*/
	void parse_input_line(std::string_view input_line)
	{
		// If line starts with #, then it's a comment and ignore
		//if (!input_line.starts_with("#"))
		if (not_comment_or_blank_line(input_line))
		{
			// We use getline except use | as our delimiter to split the string. 
			std::stringstream input_line_stream {};
			input_line_stream << input_line;
			
			std::string input_name {get_str_from_line(input_line_stream)};
			std::string input_val {get_str_from_line(input_line_stream)};

			// Save the input option into options array. 
			if(!save_input_opt(input_name, input_val))
			{
				std::cerr << "Error! Issue with this line in the input file:\n";
				std::cerr << " " << input_line << "\n";
			}

		}
	}

	/**
	* @brief Top-level function to read in options from input file
	*
	* This goes through all the input options and changes the values from 
	* their defaults using the value included in the input file. All the input 
	* options can be found in the options array at the top of this file.
	*
	* @param input_stream The filestream object pointing to the input file
	*/
	void load_input_opts(std::ifstream& input_stream)
	{
		// Load entire line as a string
		std::string input_line {};
		while (std::getline(input_stream, input_line))
		{
			// Nothing is returned by this. For each option it reads it'll just 
			// place it within the array Input::options. 
			//std::cout << "Parsing line: " << input_line << '\n';
			parse_input_line(input_line);
		}
	}
}

// Example usage.
/*
int main()
{
	std::ifstream input_stream {"testcase01.flan"};

	// If we could open the file, exit.
	if (!input_stream)
	{
		std::cerr << "Unable to open file.\n";
		return -1;
	}

	load_input_opts(input_stream);

	return 0;
}
*/
