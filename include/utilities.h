/**
* @file utilities.h
* @brief Header file for utilities.cpp
*/
#pragma once

#include <sstream>
#include <string>
#include <vector>

#include "flan_types.h"
#include "mpi.h"
#include "vectors.h"


namespace Utilities
{
	// Convert a string to an integer
	int str_as_int(const std::string& str);


	// Convert a string to a double
	double str_as_dbl(const std::string& str);


	// Split string at spaces and return it as a vector
	std::vector<std::string> split_str_at_spaces(std::string& str);


	// Cross product
	std::array<double, 3> cross_product(const std::array<double, 3>& a, 
		const std::array<double, 3>& b);


	// Dot product
	double dot_product(const std::array<double, 3>& a, 
		const std::array<double, 3>& b);


	// Bilinear interpolation at (x, y)
	double bilinear_interpolate(const double x0, const double y0, 
		const double x1, const double y1, const double z00, const double z10,
		const double z01, const double z11, const double x, const double y);


	// Trilinear interpolation at (x, y, z)
	double trilinear_interpolate(
		const double x0, const double y0, const double z0, 
		const double x1, const double y1, const double z1,
		const double v000, const double v100, const double v010, 
		const double v110, const double v001, const double v101, 
		const double v011, const double v111,
		const double x, const double y, const double z);


	// Get index of nearest neighboring cell center at val in cell idx.
	template <typename T>
	int get_neighbor_index(const double val, 
		const std::vector<T>& cell_centers, const int idx);


	// Create vector of N equally spaced values between a and b
	std::vector<double> linspace(double a, double b, std::size_t N);


	// Interpolate a Vector4D at (t0, x0, y0, z0)
	template <typename T>
	double interp_vec4d(const Vectors::Vector4D<T>& vec4d, 
		const std::vector<double> t, const std::vector<double> x,
		const std::vector<double> y, const std::vector<double> z,
		const double t0, const double x0, const double y0, const double z0);


	// Finds the indices of the two grid points that bracket a specified
	// coordinate value.
	std::pair<int, int> bracket_indices(const std::vector<double>& x, 
		double x0);


	// Broadcast vector from root rank to all other ranks
	template <typename T>
	void mpi_broadcast_vector(std::vector<T>& v, int root, MPI_Comm comm);


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
