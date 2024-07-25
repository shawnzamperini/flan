#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include "msgpack.hpp"


// function to reshape a 1D vector into a 4D
template <typename T>
std::vector<std::vector<std::vector<std::vector<T>>>> 
	reshape_4D(const std::vector<T> vec, int dim1, int dim2, int dim3, int dim4)
{
	std::vector<std::vector<std::vector<std::vector<T>>>> result(dim1,
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

int main()
{
	std::string fname {"d3d-167196-v6-gpu-elc_M0_26.gkyl"};
	std::cout << "Reading: " << fname << '\n';
	std::ifstream gkyl_stream {fname, std::ios::binary};

	// Read in the first 5 magic numbers (they are char codes for gkyl0).
	char c {};
	for (int i {}; i < 5; ++i)
	{
		gkyl_stream.read(&c, 1);
		std::cout << c << '\n';
	}

	// Next is the version, stored as 8 bytes.
	uint64_t ver {};
	gkyl_stream.read(reinterpret_cast<char*>(&ver), 8);
	std::cout << "ver = " << ver << '\n';

	// Next is the filetype
	uint64_t filetype {};
	gkyl_stream.read(reinterpret_cast<char*>(&filetype), 8);
	std::cout << "filetype = " << filetype << '\n';

	// Meta data size, number of bytes of meta data
	uint64_t meta_size {};
	gkyl_stream.read(reinterpret_cast<char*>(&meta_size), 8);
	std::cout << "meta_size = " << meta_size << '\n';

	// load	serialized data into buffer
	std::vector<char> buffer (meta_size);
	gkyl_stream.read(buffer.data(), meta_size);

	// create unpacked object
	msgpack::object_handle oh = msgpack::unpack(buffer.data(), buffer.size());

	// access root object
	msgpack::object obj = oh.get();

	// output the unpacked data
	std::cout << obj << '\n';

	std::cout << obj.type << '\n';
	std::cout << msgpack::type::MAP << '\n';

	// Access data based on the type of the object
	if (obj.type == msgpack::type::STR) {
		// If it's a string, access the string value
		std::string str = obj.as<std::string>();
		std::cout << "String value: " << str << std::endl;
	} else if (obj.type == msgpack::type::POSITIVE_INTEGER) {
		// If it's a positive integer, access the integer value
		uint64_t integer = obj.as<uint64_t>();
		std::cout << "Integer value: " << integer << std::endl;
	} else if (obj.type == msgpack::type::FLOAT32 || obj.type == msgpack::type::FLOAT64) {
		// If it's a float, access the float value
		double floating_point = obj.as<double>();
		std::cout << "Float value: " << floating_point << std::endl;
	} else if (obj.type == msgpack::type::ARRAY) {
		// If it's an array, access elements using array handle
		auto array = obj.via.array;
		for (auto const& element : array) {
			std::cout << "Element: " << element << std::endl;
		}
	} else if (obj.type == msgpack::type::MAP) {
		// If it's a map, access elements using map handle
		std::map<std::string, msgpack::type::variant> map = obj.convert();
		//auto map = obj.via.map;
		for (auto const& kv : map) {
			//std::cout << "Key: " << kv.first << ", Value: " << kv.second << std::endl;
			if (kv.second.is_string())
			{
				std::cout << "Key: " << kv.first << ", Value: " << kv.second.as_string() << std::endl;
			}
			else if (kv.second.is_double())
			{
				std::cout << "Key: " << kv.first << ", Value: " << kv.second.as_double() << std::endl;
			}
			else if (kv.second.is_int64_t())
			{
				std::cout << "Key: " << kv.first << ", Value: " << kv.second.as_int64_t() << std::endl;
			}
			else if (kv.second.is_uint64_t())
			{
				std::cout << "Key: " << kv.first << ", Value: " << kv.second.as_uint64_t() << std::endl;
			}
			//std::cout << "Key: " << kv.first << " " << kv.second.is_string() << std::endl;
		}	
	
		//std::cout << "time: " << map["time"] << '\n';
	}

	// Okay, now the header info has been loaded. What follows next depends on the file type. 
	if (filetype == 3)
	{
		// I am mostly sure that real_type = 1 is a 4-byte float, and 2 is an 8-byte float (a double).
		uint64_t real_type {};
		gkyl_stream.read(reinterpret_cast<char*>(&real_type), 8);
		std::cout << "real_type = " << real_type << '\n';

		uint64_t ndim {};
		gkyl_stream.read(reinterpret_cast<char*>(&ndim), 8);
		std::cout << "ndim = " << ndim << '\n';

		// we want to follow the pgkyl function read_domain_t1a3_v1 then read_data_t3_v1
		// grid dimensions
		//uint64_t num_dims {};
		//gkyl_stream.read(reinterpret_cast<char*>(&num_dims), 8);
		//std::cout << "num_dims = " << num_dims << '\n';

		// grid shape. this read in num_dims ints
		std::vector<uint64_t> cells {};
		uint64_t tmp_cell {};
		for (uint64_t i {}; i < ndim; ++i)
		{
			gkyl_stream.read(reinterpret_cast<char*>(&tmp_cell), 8);
			std::cout << "tmp_cell = " << tmp_cell << '\n';
			cells.push_back(tmp_cell);
		}

		// lower bounds of grid float64[ndim]
		std::cout << "sizeof(double) = " << sizeof(double) << '\n';
		double lower {};
		for (uint64_t i {}; i < ndim; ++i)
		{
			gkyl_stream.read(reinterpret_cast<char*>(&lower), 8);
			std::cout << "lower = " << lower << '\n';
		}

		// upper bounds of grid float64[ndim]
		double upper {};
		for (uint64_t i {}; i < ndim; ++i)
		{
			gkyl_stream.read(reinterpret_cast<char*>(&upper), 8);
			std::cout << "upper = " << upper << '\n';
		}

		// element size * number of components in field. divide
		// by sizeof(real_type), probably a 8 byte float, to get
		// size of elements
		uint64_t esznc {};
		gkyl_stream.read(reinterpret_cast<char*>(&esznc), 8);
		std::cout << "esznc = " << esznc << '\n';

		// number of elements is esznc / sizeof(double), so divide by 8.
		uint64_t num_comps {esznc / 8};
		std::cout << "num_comps = " << num_comps << '\n';

		// total number of cells in field
		uint64_t size {};
		gkyl_stream.read(reinterpret_cast<char*>(&size), 8);
		std::cout << "size = " << size << '\n';

		// number of ranges stored in file
		uint64_t nrange {};
		gkyl_stream.read(reinterpret_cast<char*>(&nrange), 8);
		std::cout << "nrange = " << nrange << '\n';

		// create an array for data of the final shape, gshape
		std::vector<uint64_t> gshape {};
		for (uint64_t i {}; i < ndim; ++i)
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

		// assemble the (probably 4D) array
		std::cout << "creating data vector\n";
		std::vector<std::vector<std::vector<std::vector<double>>>> data {};

		// Resize the 4D vector to match the dimensions
		data.resize(gshape[0]);
		for (uint64_t i = 0; i < gshape[0]; ++i) {
			data[i].resize(gshape[1]);
			for (uint64_t j = 0; j < gshape[1]; ++j) {
				data[i][j].resize(gshape[2]);
				for (uint64_t k = 0; k < gshape[2]; ++k) {
					data[i][j][k].resize(gshape[3]);
				}
			}
		}

		// Get the size of each dimension
		int size1 = data.size();                // Size of the 1st dimension (dim1)
		int size2 = data[0].size();             // Size of the 2nd dimension (dim2)
		int size3 = data[0][0].size();          // Size of the 3rd dimension (dim3)
		int size4 = data[0][0][0].size();       // Size of the 4th dimension (dim4)

		// Output the sizes
		std::cout << "Size of the 4D vector:" << std::endl;
		std::cout << "Dimension 1: " << size1 << std::endl;
		std::cout << "Dimension 2: " << size2 << std::endl;
		std::cout << "Dimension 3: " << size3 << std::endl;
		std::cout << "Dimension 4: " << size4 << std::endl;

		uint64_t tmp_loidx {};
		uint64_t tmp_upidx {};
		uint64_t asize {};
		std::vector<double> raw_data {};
		for (uint64_t i {}; i < nrange; ++i)
		{
			// each range starts with an uint64_t loidx and upidx which
			// are the index of lower-left and upper-right corner of the
			// range (what does this mean...?)
			std::vector<uint64_t> loidx {};
			std::vector<uint64_t> upidx {};
			for (uint64_t j {}; j < ndim; ++j)
			{
				gkyl_stream.read(reinterpret_cast<char*>(&tmp_loidx), 8);
				std::cout << "tmp_loidx = " << tmp_loidx << '\n';
				loidx.push_back(tmp_loidx);
			}
			for (uint64_t j {}; j < ndim; ++j)
			{
				gkyl_stream.read(reinterpret_cast<char*>(&tmp_upidx), 8);
				std::cout << "tmp_upidx = " << tmp_upidx << '\n';
				upidx.push_back(tmp_upidx);
			}

			// Total number of cells in range
			gkyl_stream.read(reinterpret_cast<char*>(&asize), 8);
			std::cout << "asize = " << asize << '\n';

			// read asize*esznc bytes of data
			raw_data.resize(asize * num_comps);
			for (uint64_t j {}; j < asize * num_comps; ++j)
			{
				gkyl_stream.read(reinterpret_cast<char*>(&raw_data[j]), 8);
			}

			// redefine gshape
			for (uint64_t j {}; j < ndim; ++j)
			{
				gshape[j] = upidx[j] - loidx[j] + 1;
				std::cout << "gshape[" << j << "] = " << gshape[j] << '\n';
			}

			// reshape raw_data to be put into data. after this they have the same shape
			std::cout << "raw_data.size() = " << raw_data.size() << '\n';
			auto raw_data_4d {reshape_4D(raw_data, gshape[0], gshape[1], gshape[2], gshape[3])};

			// For 4D data, loidx[d] and upidx[d] are the ranges where raw_data gets put
			// into data. This can probably be done more efficiently, but here's a simple
			// strightforward approach: loop through the whole data array, placing the
			// corresponding raw_data element into data if the index is within
			// the corresponding bounds.
			for (uint64_t i0 {loidx[0]-1}; i0 < upidx[0]; ++i0)
			{
				for (uint64_t i1 {loidx[1]-1}; i1 < upidx[1]; ++i1)
				{
					for (uint64_t i2 {loidx[2]-1}; i2 < loidx[2]; ++i2)
					{

						// We've index the part of data that this raw_data corresponds to,
						// so now we are placing raw into data.
						for (uint64_t j0 {}; j0 < raw_data_4d.size(); ++j0)
						{
							for (uint64_t j1 {}; j1 < raw_data_4d[0].size(); ++j1)
							{
								for (uint64_t j2 {}; j2 < raw_data_4d[0][0].size(); ++j2)
								{

									// Put all the components into the data array.
									//std::cout << "i0,i1,i2 = " << i0 << ", " << i1 << ", " << i2 << '\n';
									for (uint64_t c {}; c < num_comps; ++c)
									{
										data[i0][i1][i2][c] = raw_data_4d[j0][j1][j2][c];
									}
								}
							}
						}
					}
				}
			}
		}

		// At this point we have our data array (for this time slice)! Let's save it to a
		// binary file so we can validate it in python.
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
