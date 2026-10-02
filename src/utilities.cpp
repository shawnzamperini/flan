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

#include "flan_types.h"
#include "mpi.h"
#include "utilities.h"
#include "vectors.h"

namespace Utilities
{

	/**
	* Converts a string to an integer.
	*
	* The conversion is performed using std::stoi. The input string
	* must contain a valid integer representation. If the string
	* cannot be converted or the converted value falls outside the
	* range of int, std::stoi will throw an exception.
	*
	* @param str String containing the integer representation.
	* @return Integer value represented by the input string.
	* @throws std::invalid_argument If no valid conversion to int can
	*         be performed.
	* @throws std::out_of_range If the converted value is outside the
	*         range representable by int.
	*/
	int str_as_int(const std::string& str)
	{
		// Would like to add some exception handling here.
		return std::stoi(str);
	}


	/**
	* Converts a string to a double-precision floating-point value.
	*
	* The conversion is performed using std::stod. The input string
	* must contain a valid floating-point representation. If the string
	* cannot be converted or the converted value falls outside the
	* range representable by double, std::stod will throw an exception.
	*
	* @param str String containing the floating-point representation.
	* @return Double value represented by the input string.
	* @throws std::invalid_argument If no valid conversion to double can
	*         be performed.
	* @throws std::out_of_range If the converted value is outside the
	*         range representable by double.
	*/
	double str_as_dbl(const std::string& str)
	{
		// Would like to add some exception handling here.
		//std::cout << "str_as_dbl: str = " << str << '\n';
		return std::stod(str);
	}


	/**
	* Splits a string into whitespace-delimited tokens.
	*
	* Consecutive whitespace characters are treated as a single
	* delimiter, and leading or trailing whitespace is ignored.
	* The resulting vector contains each extracted token in the
	* order it appears in the input string.
	*
	* Internally, the string is parsed using a std::istringstream
	* and repeated applications of the stream extraction operator.
	*
	* @param str String to split.
	* @return Vector containing the whitespace-delimited tokens
	*         extracted from the input string.
	*/
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


	/**
	* Calculates the cross product of two three-dimensional vectors.
	*
	* Given vectors a and b, the returned vector is equal to a × b and is
	* perpendicular to both input vectors according to the right-hand rule.
	*
	* @param a First input vector.
	* @param b Second input vector.
	* @return Cross product a × b.
	*/
	std::array<double, 3> cross_product(const std::array<double, 3>& a, 
		const std::array<double, 3>& b)
	{
		return 
		{ 
			a[1] * b[2] - a[2] * b[1],
			a[2] * b[0] - a[0] * b[2],
			a[0] * b[1] - a[1] * b[0] 
		};
	}


	/**
	* Calculates the dot product of two three-dimensional vectors.
	*
	* The dot product is equal to the sum of the products of the corresponding
	* vector components.
	*
	* @param a First input vector.
	* @param b Second input vector.
	* @return Dot product a · b.
	*/
	double dot_product(const std::array<double, 3>& a, 
		const std::array<double, 3>& b)
	{
		return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	}


	/**
	* Performs bilinear interpolation within a rectangular region in
	* two-dimensional (x, y) space.
	*
	* The function interpolates a value at the coordinates (x, y)
	* using values defined at the corners of a rectangle bounded by
	* (x0, y0) and (x1, y1).
	*
	* @param x0 Lower x-coordinate of the rectangle.
	* @param y0 Lower y-coordinate of the rectangle.
	* @param z0 Value at the lower corner.
	* @param x1 Upper x-coordinate of the rectangle.
	* @param y1 Upper y-coordinate of the rectangle.
	* @param z1 Value at the upper corner.
	* @param x X-coordinate at which to evaluate the interpolated value.
	* @param y Y-coordinate at which to evaluate the interpolated value.
	* @return Bilinearly interpolated value at (x, y).
	*/
	/*
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
	*/


	/**
	 * Performs bilinear interpolation within a rectangular region in
	 * two-dimensional (x, y) space.
	 *
	 * The rectangle is defined by the coordinates (x0, y0) and
	 * (x1, y1), with values specified at each of the four corners.
	 * Interpolation is performed first along the x direction and then
	 * along the y direction to obtain the value at (x, y).
	 *
	 * Corner indexing follows the usual convention:
	 *
	 *   z00 : (x0, y0)
	 *   z10 : (x1, y0)
	 *   z01 : (x0, y1)
	 *   z11 : (x1, y1)
	 *
	 * @param x0 Lower x-coordinate of the rectangle.
	 * @param y0 Lower y-coordinate of the rectangle.
	 * @param x1 Upper x-coordinate of the rectangle.
	 * @param y1 Upper y-coordinate of the rectangle.
	 * @param z00 Value at (x0, y0).
	 * @param z10 Value at (x1, y0).
	 * @param z01 Value at (x0, y1).
	 * @param z11 Value at (x1, y1).
	 * @param x X-coordinate at which to evaluate the interpolated value.
	 * @param y Y-coordinate at which to evaluate the interpolated value.
	 * @return Bilinearly interpolated value at (x, y).
	 */
	double bilinear_interpolate(const double x0, const double y0, 
		const double x1, const double y1, const double z00, const double z10,
		const double z01, const double z11, const double x, const double y)
	{
		// Normalized coordinates.
		const double tx = (x - x0) / (x1 - x0);
		const double ty = (y - y0) / (y1 - y0);

		// Interpolate along x at y0.
		const double z_y0 = z00 + tx * (z10 - z00);

		// Interpolate along x at y1.
		const double z_y1 = z01 + tx * (z11 - z01);

		// Interpolate along y and return.
		return z_y0 + ty * (z_y1 - z_y0);
	}


	/**
	* Performs trilinear interpolation within a rectangular cell in
	* (x, y, z) space.
	*
	* The cell is defined by the coordinates (x0, y0, z0) and
	* (x1, y1, z1), with values specified at each of the eight cell
	* corners. Interpolation weights are clamped to the range [0, 1]
	* to prevent extrapolation beyond the cell boundaries. As a result,
	* points outside the cell return the value on the nearest cell face,
	* edge, or corner.
	*
	* Corner indexing follows the usual convention:
	*
	*   v000 : (x0, y0, z0)
	*   v100 : (x1, y0, z0)
	*   v010 : (x0, y1, z0)
	*   v110 : (x1, y1, z0)
	*   v001 : (x0, y0, z1)
	*   v101 : (x1, y0, z1)
	*   v011 : (x0, y1, z1)
	*   v111 : (x1, y1, z1)
	*
	* @param x0 Lower x-coordinate of the cell.
	* @param y0 Lower y-coordinate of the cell.
	* @param z0 Lower z-coordinate of the cell.
	* @param x1 Upper x-coordinate of the cell.
	* @param y1 Upper y-coordinate of the cell.
	* @param z1 Upper z-coordinate of the cell.
	* @param v000 Field value at (x0, y0, z0).
	* @param v100 Field value at (x1, y0, z0).
	* @param v010 Field value at (x0, y1, z0).
	* @param v110 Field value at (x1, y1, z0).
	* @param v001 Field value at (x0, y0, z1).
	* @param v101 Field value at (x1, y0, z1).
	* @param v011 Field value at (x0, y1, z1).
	* @param v111 Field value at (x1, y1, z1).
	* @param x X-coordinate at which to evaluate the field.
	* @param y Y-coordinate at which to evaluate the field.
	* @param z Z-coordinate at which to evaluate the field.
	* @return Interpolated field value.
	*/
	double trilinear_interpolate(
		const double x0, const double y0, const double z0, 
		const double x1, const double y1, const double z1,
		const double v000, const double v100, const double v010, 
		const double v110, const double v001, const double v101, 
		const double v011, const double v111,
		const double x, const double y, const double z)
	{
		// Normalized coordinates in [0,1]
		double tx = (x - x0) / (x1 - x0);
		double ty = (y - y0) / (y1 - y0);
		double tz = (z - z0) / (z1 - z0);

		// Prevent extrapolating past the cell center values on the edges of 
		// the grid. If we don't do this, we can get some pretty incorrect 
		// values.
		tx = std::clamp(tx, 0.0, 1.0);
		ty = std::clamp(ty, 0.0, 1.0);
		tz = std::clamp(tz, 0.0, 1.0);

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


	/**
	* Returns the index of the nearest neighboring cell center to a
	* specified cell center index.
	*
	* The neighbor is selected based on which side of the cell center
	* the value lies. If the value is greater than the cell center,
	* the neighboring index to the right is returned. If the value is
	* less than or equal to the cell center, the neighboring index to
	* the left is returned.
	*
	* At the edges of the grid, where only one neighboring cell exists,
	* that neighbor is returned regardless of which side of the cell
	* center the value lies on.
	*
	* @tparam T Data type of the cell center coordinates.
	* @param val Coordinate value of interest.
	* @param cell_centers Array of cell center coordinates.
	* @param idx Index of the reference cell center.
	* @return Index of the neighboring cell center closest to the side
	*         of the cell containing val.
	*/
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


	/**
	* Creates a vector of N equally spaced values spanning the interval
	* [a, b].
	*
	* The returned vector includes both endpoints when N > 1. If N is
	* zero, an empty vector is returned. If N is one, the returned
	* vector contains only a.
	*
	* @param a Lower bound of the interval.
	* @param b Upper bound of the interval.
	* @param N Number of values to generate.
	* @return Vector containing N equally spaced values between a and b.
	*/
	std::vector<double> linspace(double a, double b, std::size_t N)
	{
		std::vector<double> v;
		v.reserve(N);

		if (N == 0) return v;
		if (N == 1) {
			v.push_back(a);
			return v;
		}

		double step = (b - a) / (N - 1);

		for (std::size_t i = 0; i < N; ++i)
			v.push_back(a + i * step);

		return v;
	}


	/**
	* Finds the indices of the two grid points that bracket a specified
	* coordinate value.
	*
	* For values within the coordinate range, the returned indices
	* correspond to the nearest lower and upper grid points that
	* surround x0. For values outside the coordinate range, the first
	* or last valid interval is returned so that interpolation can
	* proceed without accessing elements beyond the array bounds.
	*
	* The input coordinate array is assumed to be sorted in ascending
	* order.
	*
	* @param x Sorted coordinate array.
	* @param x0 Coordinate value for which bracketing indices are sought.
	* @return Pair containing the lower and upper bracketing indices
	*         {lo, hi}.
	*/
	std::pair<int, int> bracket_indices(const std::vector<double>& x, double x0)
	{

		auto it = std::lower_bound(x.begin(), x.end(), x0);
		int i = it - x.begin();
		int N = x.size();

		// hi_raw = clamp(i, 0, N-1)
		int hi_raw = std::min(i, N - 1);

		// lo_raw = clamp(hi_raw - 1, 0, N-2)
		int lo_raw = std::max(0, hi_raw - 1);

		// mask = 1 if i == 0, else 0
		int mask = (i == 0);

		// Blend:
		//   if i == 0 → lo=0, hi=1
		//   else      → lo=lo_raw, hi=hi_raw
		int lo = lo_raw * (1 - mask) + 0 * mask;
		int hi = hi_raw * (1 - mask) + 1 * mask;

		return {lo, hi};
	}


	/**
	* Interpolates a value from a 4D field defined on a structured
	* (t, x, y, z) grid.
	*
	* Trilinear interpolation is performed in the spatial dimensions
	* (x, y, z) at the two time slices that bracket t0, followed by
	* linear interpolation in time. Spatial coordinates are clamped
	* to the grid boundaries to prevent extrapolation beyond the edge
	* cell-center values.
	*
	* @tparam T Data type stored in the 4D field.
	* @param vec4d 4D field to sample.
	* @param t Time coordinate array.
	* @param x X coordinate array.
	* @param y Y coordinate array.
	* @param z Z coordinate array.
	* @param t0 Time coordinate at which to evaluate the field.
	* @param x0 X coordinate at which to evaluate the field.
	* @param y0 Y coordinate at which to evaluate the field.
	* @param z0 Z coordinate at which to evaluate the field.
	* @return Interpolated field value.
	*/
	template <typename T>
	double interp_vec4d(const Vectors::Vector4D<T>& vec4d, 
		const std::vector<double> t, const std::vector<double> x,
		const std::vector<double> y, const std::vector<double> z,
		const double t0, const double x0, const double y0, const double z0)
	{
	
		// Similar to trilinear_interpolate, don't attempt to interpolate past
		// the last time value (since there's no way we would know what the
		// value should be). 
		const double tc = std::clamp(t0, t.front(), t.back());

		// First find the bracketing indices for each dimension
		auto [it0, it1] = bracket_indices(t, tc);
		auto [ix0, ix1] = bracket_indices(x, x0);
		auto [iy0, iy1] = bracket_indices(y, y0);
		auto [iz0, iz1] = bracket_indices(z, z0);

		// Then perform trilinear interpolation in x,y,z dimensions at each
		// time location. trilinear_interpolate prevents from interpolating
		// outside  x[ix0] to x[ix1], which can happen when you're past the
		// last cell center on the edge of the grid.
		double interp_val0 {trilinear_interpolate(
			x[ix0], y[iy0], z[iz0], 
			x[ix1], y[iy1], z[iz1], 
			vec4d(it0, ix0, iy0, iz0), vec4d(it0, ix1, iy0, iz0), // v000, v100
			vec4d(it0, ix0, iy1, iz0), vec4d(it0, ix1, iy1, iz0), // v010, v110
			vec4d(it0, ix0, iy0, iz1), vec4d(it0, ix1, iy0, iz1), // v001, v101
			vec4d(it0, ix0, iy1, iz1), vec4d(it0, ix1, iy1, iz1), // v011, v111
			x0, y0, z0)};
		double interp_val1 {trilinear_interpolate(
			x[ix0], y[iy0], z[iz0], 
			x[ix1], y[iy1], z[iz1], 
			vec4d(it1, ix0, iy0, iz0), vec4d(it1, ix1, iy0, iz0), // v000, v100
			vec4d(it1, ix0, iy1, iz0), vec4d(it1, ix1, iy1, iz0), // v010, v110
			vec4d(it1, ix0, iy0, iz1), vec4d(it1, ix1, iy0, iz1), // v001, v101
			vec4d(it1, ix0, iy1, iz1), vec4d(it1, ix1, iy1, iz1), // v011, v111
			x0, y0, z0)};

		// Then linearly interpolate in the t dimension with just point-slope
		// and return the value.
		double m {(interp_val1 - interp_val0) / (t[it1] - t[it0])};
		return m * (tc - t[it1]) + interp_val1; 
	}


	/**
	* Broadcasts a vector from a root MPI rank to all other ranks in
	* the communicator.
	*
	* The root rank provides the vector contents and size. The vector
	* size is first broadcast to all ranks, allowing non-root ranks to
	* resize their local vectors before the vector data itself is
	* broadcast.
	*
	* The element type T must have a corresponding MPI datatype defined
	* through mpi_type<T>::type.
	*
	* @tparam T Data type stored in the vector.
	* @param[in,out] v Vector to broadcast. On the root rank, contains
	*                  the source data. On all other ranks, receives the
	*                  broadcasted data.
	* @param root Rank that owns the source vector.
	* @param comm MPI communicator over which the broadcast is performed.
	*/
	template <typename T> 
	void mpi_broadcast_vector(std::vector<T>& v, int root, MPI_Comm comm)
	{
		// Get rank
	    int rank {};
		MPI_Comm_rank(comm, &rank);

		// Root sets size since it has the vector on it
		int size {};
		if (rank == root)
			size = static_cast<int>(v.size());

		// Broadcast vector size
		MPI_Bcast(&size, 1, MPI_INT, root, comm);

		// Resize vector to prepare it to recieve data from root
		if (v.size() != (size_t)size)
			v.resize(size);

		// Broadcast data to other ranks. See include/flan_types.h for an
		// explanation as to what is going on with this "mpi_type" thing, if
		// you care.
		MPI_Bcast(v.data(), size, mpi_type<T>::type, root, comm);
	}
}

// Instantiate float and double templates since we separate the declaration
// and definition
template int Utilities::get_neighbor_index<float>(
    double, const std::vector<float>&, int);
template int Utilities::get_neighbor_index<double>(
    double, const std::vector<double>&, int);

template double Utilities::interp_vec4d<float>(
	const Vectors::Vector4D<float>& vec4d, 
	const std::vector<double> t, const std::vector<double> x,
	const std::vector<double> y, const std::vector<double> z,
	const double t0, const double x0, const double y0, const double z0);
template double Utilities::interp_vec4d<double>(
	const Vectors::Vector4D<double>& vec4d, 
	const std::vector<double> t, const std::vector<double> x,
	const std::vector<double> y, const std::vector<double> z,
	const double t0, const double x0, const double y0, const double z0);

template void Utilities::mpi_broadcast_vector(std::vector<int>&, int, 
	MPI_Comm);
template void Utilities::mpi_broadcast_vector(std::vector<float>&, int, 
	MPI_Comm);
template void Utilities::mpi_broadcast_vector(std::vector<double>&, int, 
	MPI_Comm);
