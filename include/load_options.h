#pragma once

/**
* @file  load_options.h
*
* @brief Header file for load_options.cpp
*/
#include "flan_types.h" // For Inputs type
#include "options.h"

namespace Options
{
	Options load_options(const Inputs& inpts);

	void load_input(Options& opts, const Inputs& inpts);

	void crosscheck_options(Options& opts);

	template <typename T>
	bool var_correct_type(auto var, const std::string_view var_name);

	template <typename T>
	void assign_option(std::function<void(T)> set_func, auto var, 
		const std::string_view var_name);
}
