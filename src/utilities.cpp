/**
* @file utilities.cpp
* @brief Contains generic utilities used in Flan
*/
#include <array>
#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <iostream>

#include "utilities.h"

namespace Utilities
{
	int str_as_int(const std::string& str)
	{
		// Would like to add some exception handling here.
		return std::stoi(str);
	}

	double str_as_dbl(const std::string& str)
	{
		// Would like to add some exception handling here.
		//std::cout << "str_as_dbl: str = " << str << '\n';
		return std::stod(str);
	}

	std::vector<std::string> split_str_at_spaces(std::string& str)
	{
		// Split the string up into a vector of strings while getting rid of
		// repeated spaces. This can be done by casting the string as a
		// stream object, and then reading it out one word at a time with 
		// operator>>, putting the result into a vector of strings, until there
		// are no more words left.
		std::istringstream str_stream {str};
		std::vector<std::string> str_vec {};
		std::string tmp_str {};
		while (str_stream >> tmp_str)
		{
			str_vec.push_back(tmp_str);
		}
		return str_vec;
		
	}

	double bilinear_interpolate(const double x0, const double y0,
		const double z0, const double x1, const double y1, const double z1,
		const double x, const double y)
	{

		// Interpolate along x at y0
		double z_x0 {z0 + (x - x0) / (x1 - x0) * (z1 - z0)};

		// Interpolate along x at y1
		double z_x1 {z0 + (x - x0) / (x1 - x0) * (z1 - z0)};

		// Interpolate along y using the results from above and return
		return z_x0 + (y - y0) / (y1 - y0) * (z_x1 - z_x0);
	}

}
