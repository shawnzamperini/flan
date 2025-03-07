#ifndef LOAD_OPTIONS_H
#define LOAD_OPTIONS_H

#include "options.h"
#include "flan_types.h" // For Inputs type

namespace Options
{

	/**
	* @brief Controlling function for moving input options (inpts) into an
	*        Options class object.
	*
	* @param inpts Object of type @ref Inputs that is assembled in a given 
	*        simulation's input file
	*
	* @return Returns an Options object with all values from inpts written
	*         into it. If a parameter is not in inpts, the defaults are used.
	*/
	Options load_options(const Inputs& inpts);

	/**
	* @brief Loads inputs from @ref Inputs inpts object into Options opts.
	* 
	* @param opts An Options object. Generally should be passed in right after
	*        constructing it.
	* @param inpts An @ref Inputs object that is assembled in a given 
	*        simulation's input file
	*/
	void load_input(Options& opts, const Inputs& inpts);

	/**
	* @brief Checks if variant (var) is of the correct type.
	*
	* @param var A variant object that contains all the types used among the
	*            various input objects. See Inputs in \ref flan_types.h.
	* @param var_name string_view of the variant name.
	*
	* @returns Returns true if the variant (var) is of type T, false if not.
	*/
	template <typename T>
	bool var_correct_type(auto var, const std::string_view var_name);

	/**
	* @brief Assign value to Options object. Overwrites defaults.
	*
	* Uses a set function (set_func) from our Options object to 
	* set the value stored in the variant (var). If the data in var is not of
	* type T, then the value is not stored and an error is printed. 
	*/
	template <typename T>
	void assign_option(std::function<void(T)> set_func, auto var, 
		const std::string_view var_name);

}

#endif
