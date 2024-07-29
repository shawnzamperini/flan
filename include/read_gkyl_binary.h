#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <tuple>
#include "msgpack.hpp"

namespace GkylBinary
{
	// Alias this long thing for an msgpack_map type
	using Msgpack_map = std::map<std::string, msgpack::type::variant>;

	// Alias for 4D vector
	template <typename T>
	using vector4d = std::vector<std::vector<std::vector<std::vector<T>>>>;

	// Reshape a 1D vector into a 4D vector
	template <typename T>
	vector4d<T> reshape_4D(const std::vector<T> vec, int dim1, int dim2, 
		int dim3, int dim4);

	// Open a Gkeyll binary file and retun the stream object
	std::ifstream open_gkyl_file(const std::string_view fname);

	// Read the first 5 "magic numbers" from a Gkeyl file.
	void read_magic_numbers(std::ifstream& stream);

	// Read an uint64_t from stream and return it as an int.
	int read_uint64_int(std::ifstream& stream);

	// Read a double from stream and return it
	double read_double(std::ifstream& stream);

	// Read in num_vals of doubles from stream and store into a vector
	void read_double_vec(std::ifstream& stream, std::vector<double>& vec, const int num_vals);

	// Read in an msgpack MAP object from stream of byte_size bytes and return it
	Msgpack_map read_msgpack_map(std::ifstream& stream, const int byte_size);

	// Get a double from an msgpack_map object by key
	double from_msgpack_map_dbl(const Msgpack_map msgpack_map, const std::string& key);

	// Get an int from an msgpack_map object by key
	int from_msgpack_map_int(const Msgpack_map msgpack_map, const std::string& key);

	// Get a string from an msgpack_map object by key
	std::string from_msgpack_map_str(const Msgpack_map msgpack_map, 
		const std::string& key);

	// Primary function. Read in a .gkyl binary file (a frame) and return a
	// tuple of the time and the 4D data within. This should only be used with
	// move semantics. 
	std::tuple<double, vector4d<double>> load_frame(std::string_view fname);
}
