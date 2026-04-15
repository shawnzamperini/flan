/**
* @file utilities.h
* @brief Header file for utilities.cpp
*/
#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>
#include <vector>
#include <sstream>

namespace Utilities
{
	/**
	* @brief Convert string to integer
	* @param str String to convert
	* @return Return string as integer value
	*/
	int str_as_int(const std::string& str);

	/**
	* @brief Convert string to double
	* @param str String to convert
	* @return Return string as double value
	*/
	double str_as_dbl(const std::string& str);

	/**
	* @brief Convert a string into a vector of strings, splitting the string
	* apart at spaces
	* @param str The string contianing words separated by spaces
	* @return Returns vector of strings of each word
	*/
	std::vector<std::string> split_str_at_spaces(std::string& str);

	/**
	* @brief Calculate cross product
	*/
	std::array<double, 3> cross_product(const std::array<double, 3>& a, 
		const std::array<double, 3>& b);

	/**
	* @brief Calculate cross product
	*/
	double dot_product(const std::array<double, 3>& a, 
		const std::array<double, 3>& b);

	/**
	* @brief Perform a bilinear interpolation between z0=f(x0,y0) and 
	* z1=f(x1,y1) to estimate z=f(x,y).
	*
	* @param x0 Value in f(x0,y0)
	* @param y0 Value in f(x0,y0)
	* @param z0 Value at f(x0,y0)
	* @param x1 Value in f(x1,y1)
	* @param y1 Value in f(x1,y1)
	* @param z1 Value at f(x1,y1)
	* @param x Value at f(x,y)
	* @param y Value at f(x,y)
	*
	* @return Returns linearly interpolated estimate of value at f(x,y)
	*/
	double bilinear_interpolate(const double x0, const double y0,
		const double z0, const double x1, const double y1, const double z1,
		const double x, const double y);

	/**
	* @brief Perform a trilinear interpolation on a regularly spaced grid for
	* a value defined in 3D space, z=f(x,y,z).
	*
	* @param x0 Value in f(x0,y0,z0)
	* @param y0 Value in f(x0,y0,z0)
	* @param z0 Value in f(x0,y0,z0)
	* @param x1 Value in f(x1,y1,z1)
	* @param y1 Value in f(x1,y1,z1)
	* @param z1 Value in f(x1,y1,z1)
	* @param v000 Value at vertex (0,0,0), normalized coordinates
	* @param v100 Value at vertex (1,0,0), normalized coordinates
	* @param v010 Value at vertex (0,1,0), normalized coordinates
	* @param v110 Value at vertex (1,1,0), normalized coordinates
	* @param v001 Value at vertex (0,0,1), normalized coordinates
	* @param v101 Value at vertex (1,0,1), normalized coordinates
	* @param v011 Value at vertex (0,1,1), normalized coordinates
	* @param v111 Value at vertex (1,1,1), normalized coordinates
	* @param x Value at f(x,y,z)
	* @param y Value at f(x,y,z)
	* @param z Value at f(x,y,z)
	*
	* The two corners of the bounding cell are (x0,y0,z0) and (x1,y1,z1). The 
	* values at each cell vertex are v000, v100, etc., where the numbers 
	* indicate the 8 vertices in normalized coordinates.
	*
	* @return Trilinearly interpolated value from the vertex values.
	*/
	double trilinear_interpolate(
		const double x0, const double y0, const double z0, 
		const double x1, const double y1, const double z1,
		const double v000, const double v100, const double v010, 
		const double v110, const double v001, const double v101, 
		const double v011, const double v111,
		const double x, const double y, const double z);

	/**
	*
	*/
	template <typename T>
	int get_neighbor_index(const double val, 
		const std::vector<T>& cell_centers, const int idx);

	/**
	* @brief Return value at x in set of (xarr, yarr) values using linear
	* interpolation.
	* @param xarr Array of x values
	* @param yarr Array of y values
	* @param x Value to get y value at
	*
	* This is defined in the header file here because otherwise we would need
	* to have definitions in the source file for every N that we use (which
	* is of course silly and not a real option).
	*
	* @return Returns linearly interpolated value at f(x)
	*/
	template <typename T, std::size_t N>
	T linear_interpolate(const std::array<T, N>& xarr, 
		const std::array<T, N>& yarr, T x) 
		{

		// Check that arrays are of same size
		if (xarr.size() != yarr.size()) {
			throw std::invalid_argument("xarr and yarr must be of same size");
		}

		// Loop to find the interval
		for (std::size_t i = 0; i < N - 1; ++i) {
			if ((x >= xarr[i] && x <= xarr[i + 1]) 
				|| (x <= xarr[i] && x >= xarr[i + 1])) {

				// Linear interpolation
				T t = (x - xarr[i]) / (xarr[i + 1] - xarr[i]);
				return yarr[i] + t * (yarr[i + 1] - yarr[i]);
			}
		}

		std::ostringstream err {};
		err << "x value is outside the interpolation range: " << x;
		throw std::out_of_range(err.str());
	}
}

#endif
