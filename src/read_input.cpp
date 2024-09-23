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
	// Array containing all the input options along with default values.
	std::array<OptionBase*, max_input_opts> options 
	{
		new OptionStr{"gkyl_dir",            "undefined"},
		new OptionStr{"gkyl_casename",       "undefined"},
		new OptionInt{"gkyl_frame_start",              0},
		new OptionInt{"gkyl_frame_end",                1},
		new OptionStr{"gkyl_elec_name",            "elc"},
		new OptionStr{"gkyl_ion_name",             "ion"},
		new OptionStr{"gkyl_file_type",         "binary"},
		new OptionInt{"imp_atom_num",                 74},
		new OptionDbl{"imp_mass_amu",             183.84},
		new OptionInt{"imp_init_charge",               1},
		new OptionInt{"imp_num",                       1},
		new OptionDbl{"imp_xmin",                    0.0},
		new OptionDbl{"imp_xmax",                    0.0},
		new OptionStr{"imp_zstart_opt",   "single_value"},
		new OptionDbl{"imp_zstart_val",              0.0},
		new OptionStr{"imp_collisions",            "off"},
		new OptionStr{"imp_var_reduct",            "off"},
		new OptionDbl{"imp_var_reduct_freq",         0.1},
		new OptionDbl{"imp_var_reduct_min_weight",   0.1},
		new OptionDbl{"imp_time_step",              1e-7},
		new OptionDbl{"imp_source_scale_fact",		 1.0}
	};

	// Provide a string and return the corresponding enumerator in Names
	int str_to_name(std::string_view input_str)
	{
		if (input_str == "gkyl_dir")                   return gkyl_dir;
		if (input_str == "gkyl_casename")              return gkyl_casename;
		if (input_str == "gkyl_frame_start")           return gkyl_frame_start;
		if (input_str == "gkyl_frame_end")             return gkyl_frame_end;
		if (input_str == "gkyl_elec_name")             return gkyl_elec_name;
		if (input_str == "gkyl_ion_name")              return gkyl_ion_name;
		if (input_str == "gkyl_file_type")             return gkyl_file_type;
		if (input_str == "imp_atom_num")               return imp_atom_num;
		if (input_str == "imp_mass_amu")               return imp_mass_amu;
		if (input_str == "imp_init_charge")            return imp_init_charge;
		if (input_str == "imp_num")                    return imp_num;
		if (input_str == "imp_xmin")                   return imp_xmin;
		if (input_str == "imp_xmax")                   return imp_xmax;
		if (input_str == "imp_zstart_opt")             return imp_zstart_opt;
		if (input_str == "imp_zstart_val")             return imp_zstart_val;
		if (input_str == "imp_collisions")             return imp_collisions;
		if (input_str == "imp_var_reduct")             return imp_var_reduct;
		if (input_str == "imp_var_reduct_freq")        return imp_var_reduct_freq;
		if (input_str == "imp_var_reduct_min_weight")  return imp_var_reduct_min_weight;
		if (input_str == "imp_time_step")              return imp_time_step;
		if (input_str == "imp_source_scale_fact")      return imp_source_scale_fact;
		else return -1;
	}

	// Should be as many Names as there are options in the array.
	// Note: It is on the developer that the order of the array and the names
	// in Names are the same! The actual input file can be in whatever order. 
	static_assert(std::ssize(options) == max_input_opts);

	// Helper function to provide a Name from the enumeration above and it will
	// return the value held by the option. One function for strings, integers
	// and doubles. 
	// The reason we do not use a template here is because the type to be 
	// returned can not be deduced with CTAD.
	std::string_view get_opt_str(Names name)
	{
		OptionBase* opt_base_ptr {options[name]};
		return dynamic_cast<OptionStr*>(opt_base_ptr)->get_value();
	}

	int get_opt_int(Names name)
	{
		OptionBase* opt_base_ptr {options[name]};
		return dynamic_cast<OptionInt*>(opt_base_ptr)->get_value();
	}

	double get_opt_dbl(Names name)
	{
		OptionBase* opt_base_ptr {options[name]};
		return dynamic_cast<OptionDbl*>(opt_base_ptr)->get_value();
	}

	// Takes in a line from an input file as a stringstream object and returns
	// it as a normal string for further manipulation. 
	std::string get_str_from_line(std::stringstream& input_line_stream)
	{
		std::string str {};
		std::getline(input_line_stream, str, '|');
		std::string::iterator end_pos = std::remove(str.begin(), str.end(), ' ');
		str.erase(end_pos, str.end());
		return str;
	}

	// Using the option name (input_name) and the value pulled from the input
	// file (input_val), save this value into the main options array. Returns
	// true if successful, false if not.
	bool save_input_opt(std::string_view input_name, std::string_view input_val)
	{
			// Assign each variable to its appropriate spot in options.
			//std::cout << Input::str_to_name(input_name) << '\n';
			Input::OptionBase* opt_base_ptr {Input::options[Input::str_to_name(input_name)]};

			// Attempt to downcast to each inherited class until we hit upon the
			// right one. Once we do, call put_value to overwrite the default value.
			if (dynamic_cast<Input::OptionStr*>(opt_base_ptr))
			{
				//std::cout << "String option\n";
				dynamic_cast<Input::OptionStr*>(opt_base_ptr)->put_value(input_val);
				//std::cout << *dynamic_cast<Input::OptionStr*>(opt_base_ptr) << '\n';
			}
			else if (dynamic_cast<Input::OptionInt*>(opt_base_ptr))
			{
				//std::cout << "Int option\n";
				dynamic_cast<Input::OptionInt*>(opt_base_ptr)->put_value(input_val);
				//std::cout << *dynamic_cast<Input::OptionInt*>(opt_base_ptr) << '\n';
			}
			else if (dynamic_cast<Input::OptionDbl*>(opt_base_ptr))
			{
				//std::cout << "Double option\n";
				dynamic_cast<Input::OptionDbl*>(opt_base_ptr)->put_value(input_val);
				//std::cout << *dynamic_cast<Input::OptionDbl*>(opt_base_ptr) << '\n';
			}
			else
			{
				std::cout << "Error: could not save input value\n";
				return false;
			}

			return true;
	}

	bool not_comment_or_blank_line(const std::string_view line) 
	{
		return !line.starts_with("#") && 
			(line.find_first_not_of(' ') != std::string::npos);
	}

	// Take a line from an input file and parse it, saving the input into 
	// the options array if the line is not blank or a comment. 
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

			// First call grabs the name, second call the value (as a string).
			//std::getline(input_line_stream, input_name, '|');
			//std::cout << "Input name: " << input_name << '\n';
			//std::getline(input_line_stream, input_val, '|');
			//std::cout << "Input value: " << input_val << '\n';

			// Save the input option into options array. 
			if(!save_input_opt(input_name, input_val))
			{
				std::cerr << "Error! Issue with this line in the input file:\n";
				std::cerr << " " << input_line << "\n";
			}

		}
	}

	// This goes through all the input options and changes the values from 
	// their defaults using the value included in the input file. All the input 
	// options can be found in the options array.
	void load_input_opts(std::ifstream& input_stream)
	{
		// Load entire line as a string
		std::string input_line {};
		while (std::getline(input_stream, input_line))
		{
			// Nothing is returned by this. For each option it reads it'll just 
			// place it within the array Input::settings. 
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
