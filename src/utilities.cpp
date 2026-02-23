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

	double trilinear_interpolate(
		const double x0, const double y0, const double z0, 
		const double x1, const double y1, const double z1,
		const double v000, const double v100, const double v010, 
		const double v110, const double v001, const double v101, 
		const double v011, const double v111,
		const double x, const double y, const double z)
	{
		// Normalized coordinates in [0,1]
		const double tx = (x - x0) / (x1 - x0);
		const double ty = (y - y0) / (y1 - y0);
		const double tz = (z - z0) / (z1 - z0);

		// Interpolate along x for the four lower/upper face corners
		const double c00 = v000 + tx * (v100 - v000);
		const double c01 = v001 + tx * (v101 - v001);
		const double c10 = v010 + tx * (v110 - v010);
		const double c11 = v011 + tx * (v111 - v011);

		// Interpolate along y for the lower and upper edges
		const double c0 = c00 + ty * (c10 - c00);
		const double c1 = c01 + ty * (c11 - c01);

		// Interpolate along z and return
		return c0 + tz * (c1 - c0);
	}

	// Find nearest neighbor index by essentially seeing which half of the
	// cell a particle is in and returning the index of the neighboring cell
	// closest to that half. Designed to be run once per x,y,z.
	template <typename T>
	int get_neighbor_index(const double val, 
		const std::vector<T>& cell_centers, const int idx)
	{
		// dx > 0 --> (1*2 - 1) = +1
		// dx < 0 --> (0*2 - 1) = -1
		double dx_from_center {val - cell_centers[idx]};
		int side {2 * (dx_from_center > 0.0) - 1};

		// Need to check we aren't at grid edges
		int is_left_edge {(idx == 0)};
		int is_right_edge {(idx == std::ssize(cell_centers) - 1)};

		// This will correctly do -1 --> 1 if we're at the left edge, and
		// 1 --> -1 if we're at the right edge. This effectively means we
		// are using the only neighboring option.
		int offset {
			  side * (1 - is_left_edge - is_right_edge)
			+ 1    * is_left_edge
			+ (-1) * is_right_edge};

		return idx + offset;
	}
}

// Instatiate float and double templates since we separate the declaration
// and definition
template int Utilities::get_neighbor_index<float>(
    double, const std::vector<float>&, int);
template int Utilities::get_neighbor_index<double>(
    double, const std::vector<double>&, int);

