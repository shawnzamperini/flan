#pragma once
/**
* @file utilities.h
* @brief Header file for utilities.cpp
*
* The linear interpolation templates for fixed-size tables are defined here,
* since they work for any table length N. The CUDA versions are in
* cuda/utilities.cuh.
*/

#include <algorithm>
#include <array>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "flan_types.h"
#include "mpi.h"
#include "vectors.h"


namespace Utilities
{
	int str_as_int(const std::string& str);

	double str_as_dbl(const std::string& str);

	std::vector<std::string> split_str_at_spaces(std::string& str);

	std::array<double, 3> cross_product(const std::array<double, 3>& a, 
		const std::array<double, 3>& b);

	double dot_product(const std::array<double, 3>& a, 
		const std::array<double, 3>& b);

	double bilinear_interpolate(const double x0, const double y0, 
		const double x1, const double y1, const double z00, const double z10,
		const double z01, const double z11, const double x, const double y);

	double trilinear_interpolate(
		const double x0, const double y0, const double z0, 
		const double x1, const double y1, const double z1,
		const double v000, const double v100, const double v010, 
		const double v110, const double v001, const double v101, 
		const double v011, const double v111,
		const double x, const double y, const double z);

	template <typename T>
	int get_neighbor_index(const double val, 
		const std::vector<T>& cell_centers, const int idx);

	std::vector<double> linspace(double a, double b, std::size_t N);

	template <typename T>
	double interp_vec4d(const Vectors::Vector4D<T>& vec4d, 
		const std::vector<double> t, const std::vector<double> x,
		const std::vector<double> y, const std::vector<double> z,
		const double t0, const double x0, const double y0, const double z0);

	std::pair<int, int> bracket_indices(const std::vector<double>& x, 
		double x0);

	template <typename T>
	void mpi_broadcast_vector(std::vector<T>& v, int root, MPI_Comm comm);

	/**
	* Linearly interpolates a fixed-size table at x.
	*
	* Searches the table from the start for the interval containing x, so
	* xarr can be increasing or decreasing. For a large, evenly spaced
	* table, linear_interpolate_uniform is much faster.
	*
	* This is defined in the header file here because otherwise we would need
	* to have definitions in the source file for every N that we use (which
	* is of course silly and not a real option).
	*
	* @tparam T Floating-point type of the table.
	* @tparam N Length of the table, deduced from the arrays.
	* @param xarr Table x values.
	* @param yarr Table y values.
	* @param x Point to interpolate at.
	* @return Interpolated y value.
	* @throws std::out_of_range If x is outside the range of xarr.
	*/
	template <typename T, std::size_t N>
	T linear_interpolate(const std::array<T, N>& xarr, 
		const std::array<T, N>& yarr, T x) 
		{

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


	/**
	* Linearly interpolates a fixed-size table with evenly spaced x values
	* at x.
	*
	* Faster than linear_interpolate for large tables, since the interval
	* containing x is calculated directly instead of searched for. xarr must
	* be increasing and evenly spaced, like the Nanbu table from
	* python/calc_nanbu_A.py. Points outside the table are extrapolated from
	* the first or last interval. Same as linear_interpolate_uniform_cuda
	* in cuda/utilities.cuh.
	*
	* @tparam T Floating-point type of the table.
	* @tparam N Length of the table, deduced from the arrays.
	* @param xarr Table x values, increasing and evenly spaced.
	* @param yarr Table y values.
	* @param x Point to interpolate at.
	* @return Interpolated y value.
	*/
	template <typename T, std::size_t N>
	T linear_interpolate_uniform(const std::array<T, N>& xarr,
		const std::array<T, N>& yarr, T x)
	{
		// Table spacing, from the ends so rounding in the individual
		// entries averages out
		const T dx {(xarr[N - 1] - xarr[0]) / (N - 1)};

		// Interval containing x. Clamped so x at (or rounding past) the
		// ends of the table stays inside it.
		int i {static_cast<int>((x - xarr[0]) / dx)};
		i = std::clamp(i, 0, static_cast<int>(N) - 2);

		// Use the actual interval width, since the table entries are
		// rounded and not exactly dx apart
		const T t {(x - xarr[i]) / (xarr[i + 1] - xarr[i])};
		return yarr[i] + t * (yarr[i + 1] - yarr[i]);
	}

}  // namespace Utilities
