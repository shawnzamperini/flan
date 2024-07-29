#ifndef READ_INPUT_H
#define READ_INPUT_H

#include <array>
#include "input_classes.h"

/*
This contains the declarations of everything needed to read in input from
an input file into Option classes used by Flan to keep track of things.
An input file takes the form of, e.g.,

# This is a comment
string_option_name  |  name     Note: Do not use quotation marks.
int_option_name     |  3
double_option_name  |  183.34

The options can go in any order, and you can move the | around to make 
everything align if you want to. The options must be valid options though or
it will yell at you!

For a list of all the options and their defaults, see the options array in
read_input.cpp (or the documentation). 
*/
namespace Input
{
	
	// Enumeration containing all the input options. There are used to index
	// the options array defined in read_input.cpp. 
	enum Names
	{
		gkyl_casename,
		gkyl_frame_start,
		gkyl_frame_end,
		imp_mass,
		max_input_opts,
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
	// file (input_val), save this value into the main options array. 
	void save_input_opt(std::string_view input_name, std::string_view input_val);

	// Take a line from an input file and parse it, saving the input into 
	// the options array if the line is not blank or a comment. 
	void parse_input_line(std::string_view input_line);

	// This goes through all the input options and changes the values from 
	// their defaults using the value included in the input file. All the input 
	// options can be found in the options array.
	void load_input_opts(std::ifstream& input_stream);
}

#endif
