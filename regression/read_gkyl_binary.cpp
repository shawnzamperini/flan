#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <type_traits>
#include <variant>
#include "msgpack.hpp"

template <typename T>
using vector4d = std::vector<std::vector<std::vector<std::vector<T>>>>;

// function to reshape a 1D vector into a 4D
template <typename T>
vector4d<T> reshape_4D(const std::vector<T> vec, int dim1, int dim2, int dim3, 
					   int dim4)
{
	vector4d result(dim1,
		std::vector<std::vector<std::vector<T>>>(dim2,
        std::vector<std::vector<T>>(dim3,
		std::vector<T>(dim4))));

	// Copy elements from the original vector to the reshaped 4D vector
    int index = 0;
    for (int i = 0; i < dim1; ++i) {
        for (int j = 0; j < dim2; ++j) {
            for (int k = 0; k < dim3; ++k) {
                for (int l = 0; l < dim4; ++l) {
                    result[i][j][k][l] = vec[index++];
                }
            }
        }
    }

    return result;
}

std::ifstream open_gkyl_file(const std::string& fname)
{
	//std::string fname {"d3d-167196-v6-gpu-elc_M0_26.gkyl"};
	std::cout << "Reading: " << fname << '\n';
	std::ifstream stream {fname, std::ios::binary};

	if (!stream.is_open())
	{
		std::cerr << "Error! Could not open file: " << fname << '\n';
	}

	return stream;
}

void read_magic_numbers(std::ifstream& stream)
{
	// The first 5 chars are "magic numbers" (they are char codes for the 
	// word "gkyl0").
	char c {};
	std::string magic_in {};
	std::string magic {"gkyl0"};
	for (int i {}; i < 5; ++i)
	{
		stream.read(&c, 1);
		magic_in.push_back(c);
	}

	// Check that the correct chars (or string) was read in.
	if (magic_in != magic)
	{
		std::cerr << "Error! Magic characters read in do not match " 
			"\"gkyl0\": " << magic_in << '\n';
	}
}

// Read an unsigned 64-bit int from stream and return as an int
int read_uint64_int(std::ifstream& stream)
{
	uint64_t var {};

	// Read in 8 bytes (64 bits) into an uint64_t
	stream.read(reinterpret_cast<char*>(&var), 8);
	return static_cast<int>(var);
}

// Read a double from stream
double read_double(std::ifstream& stream)
{
	double var {};

	// Read in 8 bytes (64 bits) into an double
	stream.read(reinterpret_cast<char*>(&var), 8);
	return var;
}

void read_double_vec(std::ifstream& stream, std::vector<double>& vec, const int num_vals)
{
	// Ensure vector is the right size (this is probably unnecessary for how 
	// we use this, but the performance penalty is negligible compared to
	// the safety it brings).
	vec.resize(num_vals);
	for (int i {}; i < num_vals; ++i)
	{
		vec[i] = read_double(stream);
	}
}

// Alias this long thing for an msgpack_map type
typedef std::map<std::string, msgpack::type::variant> Msgpack_map;

Msgpack_map read_msgpack_map(std::ifstream& stream, const int byte_size)
{
	// Load	serialized data into buffer
	std::vector<char> buffer (byte_size);
	stream.read(buffer.data(), byte_size);

	// Create unpacked object
	msgpack::object_handle oh = msgpack::unpack(buffer.data(), buffer.size());

	// Access root object
	msgpack::object obj = oh.get();

	// Output the unpacked data (debugging)
	std::cout << obj << '\n';
	std::cout << obj.type << '\n';
	std::cout << msgpack::type::MAP << '\n';

	// Ensure the data is a map type. SHould add error detection here.
	if (obj.type != msgpack::type::MAP) 
	{
		std::cerr << "Error! msgpack header info is not a map type.\n";
		std::cerr << "  obj.type = " << obj.type << '\n';
		std::cerr << "  msgpack::type::MAP" << msgpack::type::MAP << '\n';
	}	
	
	return obj.convert();
}

// Return a double entry from an msgpack_map object
double from_msgpack_map_dbl(const Msgpack_map msgpack_map, const std::string& key)
{
	for (auto const& kv : msgpack_map) 
	{
		if (kv.first == key)
		{
			if (kv.second.is_double())
			{
				std::cout << "Key: " << kv.first << ", Value: " << 
					kv.second.as_double() << std::endl;
				return kv.second.as_double();
			}
			else
			{
				std::cerr << "Error! " << key << " is not a double in "
					"msgpack_map.\n";
				return 0.0;
			}
		}	
	}

	std::cerr << "Error! " << key << " was not found in msgpack_map.\n";
	return 0.0;
}

// Return an integer entry from an msgpack_map object
int from_msgpack_map_int(const Msgpack_map msgpack_map, const std::string& key)
{
	for (auto const& kv : msgpack_map) 
	{
		if (kv.first == key)
		{
			if (kv.second.is_uint64_t())
			{
				std::cout << "Key: " << kv.first << ", Value: " << 
					kv.second.as_uint64_t() << std::endl;
				return static_cast<int>(kv.second.as_uint64_t());
			}
			else if (kv.second.is_int64_t())
			{
				std::cout << "Key: " << kv.first << ", Value: " << 
					kv.second.as_int64_t() << std::endl;
				return static_cast<int>(kv.second.as_int64_t());
			}
			else
			{
				std::cerr << "Error! " << key << " is not a int64_t/uint64_t "
					"in msgpack_map.\n";
				return 0;
			}
		}	
	}

	std::cerr << "Error! " << key << " was not found in msgpack_map.\n";
	return 0;
}

// Return a string entry from an msgpack_map object
std::string from_msgpack_map_str(const Msgpack_map msgpack_map, 
	const std::string& key)
{
	for (auto const& kv : msgpack_map) 
	{
		if (kv.first == key)
		{
			if (kv.second.is_string())
			{
				std::cout << "Key: " << kv.first << ", Value: " << 
					kv.second.as_string() << std::endl;
				return kv.second.as_string();
			}
			else
			{
				std::cerr << "Error! " << key << " is not a string in "
					"msgpack_map.\n";
				return "";
			}
		}	
	}

	std::cerr << "Error! " << key << " was not found in msgpack_map.\n";
	return "";
}


int main()
{
	// We make the assumption that a double is 8 bytes just to not get bogged
	// down in implementation details in the start. If you encounter this
	// error, it's easily fixed just ask Shawn.
	if (sizeof(double) != 8)
	{
		std::cerr << "Error! double is not 8 bytes on this machine. Ask"
			" Shawn to fix this\n";
		return -1;
	}

	// Load file (stream)
	std::ifstream gkyl_stream {open_gkyl_file("d3d-167196-v6-gpu-elc_M0_26.gkyl")};

	// Read the first 5 chars, i.e., the "magic numbers"
	read_magic_numbers(gkyl_stream);

	// Next is the version
	int version {read_uint64_int(gkyl_stream)};
	std::cout << "version = " << version << '\n';

	// Next is the filetype
	int file_type {read_uint64_int(gkyl_stream)};
	std::cout << "file_type = " << file_type << '\n';

	// Meta data size, number of bytes of meta data
	int meta_size {read_uint64_int(gkyl_stream)};
	std::cout << "meta_size = " << meta_size << '\n';

	// The next section of data is stored in an msgpack map object
	Msgpack_map msgpack_map {read_msgpack_map(gkyl_stream, meta_size)};

	// Load info from the msgpack_map into normal variables
	double time {from_msgpack_map_dbl(msgpack_map, "time")};
	int frame {from_msgpack_map_int(msgpack_map, "frame")};
	int poly_order {from_msgpack_map_int(msgpack_map, "polyOrder")};
	std::string basis_type {from_msgpack_map_str(msgpack_map, "basisType")};

	// Okay, now the header info has been loaded. What follows next depends on 
	// the file type. 
	if (file_type == 3)
	{
		// real_type = 1 is a 4-byte float, and 2 is an 8-byte float (a double).
		int real_type {read_uint64_int(gkyl_stream)};
		std::cout << "real_type = " << real_type << '\n';

		int ndim {read_uint64_int(gkyl_stream)};
		std::cout << "ndim = " << ndim << '\n';

		// grid shape. this read in ndim ints
		std::vector<int> cells {};
		int tmp_cell {};
		for (int i {}; i < ndim; ++i)
		{
			tmp_cell = read_uint64_int(gkyl_stream);
			cells.push_back(tmp_cell);
		}

		// lower bounds of grid float64[ndim]
		std::vector<double> lower {};
		for (int i {}; i < ndim; ++i)
		{
			lower.push_back(read_double(gkyl_stream));
		}

		// upper bounds of grid float64[ndim]
		std::vector<double> upper {};
		for (int i {}; i < ndim; ++i)
		{
			upper.push_back(read_double(gkyl_stream));
		}

		// element size * number of components in field. divide
		// by sizeof(real_type), probably a 8 byte float, to get
		// size of elements
		int esznc {read_uint64_int(gkyl_stream)};
		std::cout << "esznc = " << esznc << '\n';

		// number of elements is esznc / sizeof(double), so divide by 8.
		int num_comps {esznc / 8};
		std::cout << "num_comps = " << num_comps << '\n';

		// total number of cells in field
		int size {read_uint64_int(gkyl_stream)};
		std::cout << "size = " << size << '\n';

		// number of ranges stored in file
		int nrange {read_uint64_int(gkyl_stream)};
		std::cout << "nrange = " << nrange << '\n';

		// gshape is going to be the shape of the data array that we eventually
		// return. it contains the number of cells in each direction and the
		// last element is the number of components.
		std::vector<int> gshape {};
		for (int i {}; i < ndim; ++i)
		{
			gshape.push_back(cells[i]);
		}
		gshape.push_back(num_comps);
		std::cout << "gshape: ";
		for (auto g : gshape)
		{
			std::cout << g << " ";
		}
		std::cout << '\n';

		// Assemble the (probably 4D) array
		std::cout << "creating data vector\n";
		vector4d<double> data {};

		// Resize the 4D vector to match the dimensions
		data.resize(gshape[0]);
		for (int i = 0; i < gshape[0]; ++i) {
			data[i].resize(gshape[1]);
			for (int j = 0; j < gshape[1]; ++j) {
				data[i][j].resize(gshape[2]);
				for (int k = 0; k < gshape[2]; ++k) {
					data[i][j][k].resize(gshape[3]);
				}
			}
		}

		// Get the size of each dimension
		int size1 = data.size();
		int size2 = data[0].size();
		int size3 = data[0][0].size();
		int size4 = data[0][0][0].size();

		// Output the sizes
		std::cout << "Size of the 4D vector:" << std::endl;
		std::cout << "Dimension 1: " << size1 << std::endl;
		std::cout << "Dimension 2: " << size2 << std::endl;
		std::cout << "Dimension 3: " << size3 << std::endl;
		std::cout << "Dimension 4: " << size4 << std::endl;

		std::vector<double> raw_data {};
		for (int i {}; i < nrange; ++i)
		{
			// each range starts with an uint64_t loidx and upidx which
			// are the index of lower-left and upper-right corner of the
			// range (what does this mean...?)
			std::vector<int> loidx {};
			std::vector<int> upidx {};
			for (int j {}; j < ndim; ++j)
			{
				loidx.push_back(read_uint64_int(gkyl_stream));
			}
			for (int j {}; j < ndim; ++j)
			{
				upidx.push_back(read_uint64_int(gkyl_stream));
			}

			// Total number of cells in range
			int asize {read_uint64_int(gkyl_stream)};
			std::cout << "asize = " << asize << '\n';

			// read asize*esznc bytes of data
			//raw_data.resize(asize * num_comps);
			//for (int j {}; j < asize * num_comps; ++j)
			//{
			//	gkyl_stream.read(reinterpret_cast<char*>(&raw_data[j]), 8);
			//}
			read_double_vec(gkyl_stream, raw_data, asize * num_comps);

			// Redefine gshape. This will be the shape of raw_data when we
			// reshape it
			for (int j {}; j < ndim; ++j)
			{
				gshape[j] = upidx[j] - loidx[j] + 1;
				std::cout << "gshape[" << j << "] = " << gshape[j] << '\n';
			}

			// reshape raw_data so it can be mapped to data
			std::cout << "raw_data.size() = " << raw_data.size() << '\n';
			vector4d<double> raw_data_4d {reshape_4D(raw_data, gshape[0], 
				gshape[1], gshape[2], gshape[3])};

			// For 4D data, loidx[d] and upidx[d] are the ranges where raw_data 
			// gets put into data. So loop through all those elements in data
			// and place the corresponding raw_data component into data
			for (int i0 {loidx[0]-1}; i0 < upidx[0]; ++i0)
			{
				for (int i1 {loidx[1]-1}; i1 < upidx[1]; ++i1)
				{
					for (int i2 {loidx[2]-1}; i2 < upidx[2]; ++i2)
					{

						// Put all the components into the data array. raw_data
						// is not the same size as data, but you can think of
						// it as a subset of data that has been reindexed to
						// start at zero in each dimension. loidx[0]-1 is just
						// how we access an element in raw_data that matches
						// data.
						//std::cout << "i0,i1,i2 = " << i0 << ", " << i1 << 
						//	", " << i2 << '\n';
						for (int c {}; c < num_comps; ++c)
						{
							data[i0][i1][i2][c] = raw_data_4d[i0-(loidx[0]-1)]
								[i1-(loidx[1]-1)][i2-(loidx[2]-1)][c];
						}
					}
				}
			}
		}

		// At this point we have our data array (for this time slice)! Let's 
		// save it to a binary file so we can validate it in python.
		std::ofstream outfile {"data.bin", std::ios::binary};

		// Save the dimensions
		std::cout << "Writing array dimensions...\n";
		uint64_t d0 {std::size(data)};
		uint64_t d1 {std::size(data[0])};
		uint64_t d2 {std::size(data[0][0])};
		uint64_t d3 {std::size(data[0][0][0])};
		outfile.write(reinterpret_cast<const char*>(&d0), sizeof(uint64_t));
		outfile.write(reinterpret_cast<const char*>(&d1), sizeof(uint64_t));
		outfile.write(reinterpret_cast<const char*>(&d2), sizeof(uint64_t));
		outfile.write(reinterpret_cast<const char*>(&d3), sizeof(uint64_t));

		// Write out array
		std::cout << "Writing array...\n";
		for (const auto& v0 : data)
		{
			for (const auto& v1 : v0)
			{
				for (const auto& v2 : v1)
				{
					outfile.write(reinterpret_cast<const char*>(v2.data()), v2.size() * sizeof(double));
				}
			}
		}

		outfile.close();

	}


	gkyl_stream.close();

	return 0;
}
