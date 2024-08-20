#include <iostream>
#include <string>
#include <sstream>
#include <tuple>
#include <type_traits> // For std::is_same
#include <numeric>     // For std::accumulate
#include <functional>  // For std::multiplies
/*
 * Opting out of ADIOS2 since I think it's going to get removed...
#include <adios2.h>
*/
#include "read_input.h"
#include "read_gkyl.h"
//#include "read_gkyl_binary.h"
#include "vectors.h"
#include "constants.h"
#include "background.h"

namespace Gkyl
{
	// Vectors to hold density, temperature, potential and magnetic field
	// for each frame (assuming electrostatic so only one magnetic field
	// entry is needed. Dimensions are (time, x, y, z). 
	Vectors::Vector4D gkyl_ne {};
	Vectors::Vector4D gkyl_te {};
	Vectors::Vector4D gkyl_ti {};
	Vectors::Vector4D gkyl_vp {};
	Vectors::Vector4D gkyl_b {};  // If electrostatic the first dimension is only 1 long

	// Entry point for reading Gkeyll data into Flan. 
	Background::Background read_gkyl()
	{	
		// Load each needed dataset from Gkeyll
		read_elec_density();
		read_elec_temperature();
		read_ion_temperature();
		read_potential();
		read_magnetic_field();

		// With all our arrays assembled, encapsulate them into a Background
		// class object with move semantics then return.
		Background::Background bkg {create_bkg()};
		return bkg;
	}
	
	// Simple helper function to return a string of the full path to a Gkeyll
	// file using it's file name convention.
	std::string assemble_path(const std::string& species, 
		const std::string& ftype, int frame, const std::string& extension) 
	{
		// Load into local variables so code is easier to read.
		std::string gkyl_dir {Input::get_opt_str(Input::gkyl_dir)};
		std::string gkyl_casename {Input::get_opt_str(Input::gkyl_casename)};

		return gkyl_dir + "/" + gkyl_casename + "-" + species + "_" + ftype + 
			"_" + std::to_string(frame) + extension;
	}

	// Function to return the full path to read_gkyl.py, which is used to
	// interface with postgkyl and create files that are easily read in by
	// Flan. This path is stored as an environment variable that is set
	// before Flan is run.
	std::string get_read_gkyl_py()
	{
		if (const char* read_gkyl_py {std::getenv("READ_GKYL_PY")})
		{
			return static_cast<std::string>(read_gkyl_py);
		}
		else
		{
			std::cerr << "Error! READ_GKYL_PY is not set. Please set this "
				<< "to the full path to the read_gkyl.py script. E.g.:\n"
				<< "  READ_GKYL_PY=/path/to/flan/python/read_gkyl.py\n";
			return std::string {};
		}
	}

	// Reads in a 1D vector of times corresponding to each frame.
	std::vector<double> load_times()
	{
		// Filename is treated as a constant.
		std::string filename {"bkg_from_pgkyl_times.csv"};
		std::ifstream fstream {filename};

		// Go through one line at a time
		std::string line {};
		std::vector<double> times {};
		int ntimes {};
		bool ntimes_read {false};
		while(std::getline(fstream, line))
		{
			// If line starts with #, then it's a comment and ignore
			if (!line.starts_with("#"))
			{
				// The first line is a single integer telling us how many
				// times (frames) there are.
				if (!ntimes_read)
				{
					ntimes = std::stoi(line);
					ntimes_read = true;
					std::cout << "Reading " << ntimes << " times\n";

					// Set the capacity of the vector now that we know what 
					// it will be.
					times.reserve(ntimes);
				}
				else
				{
					// Append time to our vector. Note we already set
					// the capacity of the vector, so this should not be too
					// expensive.  
					times.push_back(std::stod(line));	
				}
			}
		}

		// Close and return (uses move semantics)
		fstream.close();
		return times;
	}

	// Reads in the x, y, z grid nodes using pgkyl.
	std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> 
		load_grid()
	{
		// Filename is treated as a constant.
		std::string filename {"bkg_from_pgkyl_grid.csv"};
		std::ifstream fstream {filename};

		// Define variables needed for the following loop
		std::string line {};
		std::vector<double> grid_x {};
		std::vector<double> grid_y {};
		std::vector<double> grid_z {};
		std::vector<int> dims {};
		int ngrid {};
		bool ngrid_read {false};
		int count {};

		// Go through one line at a time
		while(std::getline(fstream, line))
		{
			// If line starts with #, then it's a comment and ignore
			if (!line.starts_with("#"))
			{
				// The first line is a three integers telling us how many
				// the x, y, z dimensions. These are the edges of each grid
				// cell, so there is one extra value in each dimension when
				// compared to the actual value (e.g., density) dimensions.
				if (!ngrid_read)
				{
					// Use string stream to pull the three integers out.
					std::stringstream ss {line};
					int tmp_int {};
					while (ss >> tmp_int)
					{
						dims.push_back(tmp_int);	
					}

					// If we have anything other than 3 values something
					// went wrong.
					if (dims.size() != 3)
					{
						std::cerr << "Error! More than 3 dimensions were read"
							<< " in from " << filename << ". Check that this "
							<< " file is being generated correctly.\n";
						std::cerr << "dims: ";
						for (auto i : dims) std::cerr << i << " ";
						std::cerr << '\n';
					}
				
					// ngrid is the total number of points in the file. The
					// x values are all printed, then the y then the z. 
					ngrid = dims[0] + dims[1] + dims[2];
					ngrid_read = true;
					std::cout << "Reading " << ngrid << " grid values\n";

					// Set the capacity of the vectors now that we know them.
					grid_x.reserve(dims[0]);
					grid_y.reserve(dims[1]);
					grid_z.reserve(dims[2]);
				}
				else
				{
					// Append to each vector accordingly
					double grid_value {std::stod(line)};
					if (count < dims[0]) 
						grid_x.push_back(grid_value);
					else if (count >= dims[0] && count < dims[0] + dims[1]) 
						grid_y.push_back(grid_value);
					else if (count >= dims[0] + dims[1]) 
						grid_z.push_back(grid_value);

					// If count ends up being greater than all the values
					// we need to read in, something went wrong.
					if (count >= dims[0] + dims[1] + dims[2])
					{
						std::cerr << "Error! There appears to be too many "
							<< "values in " << filename << ". Check that "
							<< "this file is being generated correctly.\n";
						std::cerr << "count = " << count << "dims = "
							<< dims[0] << ", " << dims[1] << ", " << dims[2]
							<< '\n';
					}

					++count;
				}
			}
		}
		fstream.close();

		// Pack up into tuple and return
		return std::make_tuple(grid_x, grid_y, grid_z);
	}

	// Read in data values using pgkyl, returning as a Vector4D.
	Vectors::Vector4D load_values(const std::string& data_type)
	{
		// Filename is treated as a constant.
		std::string filename {"bkg_from_pgkyl_" + data_type + ".csv"};
		std::ifstream fstream {filename};

		// Define variables needed for the following loop
		std::string line {};
		std::vector<double> data_flattened {};
		std::vector<int> dims {};
		int ndata {};
		bool ndata_read {false};
		int count {};

		// Go through one line at a time
		while(std::getline(fstream, line))
		{
			// If line starts with #, then it's a comment and ignore
			if (!line.starts_with("#"))
			{
				// The first number is how many times (frames) there
				// are then the x, y, z sizes. 
				// The x, y, z dimensions are one smaller than grid_x, y, z
				// because these are the data values at the cell centers.
				if (!ndata_read)
				{
					// Use string stream to pull the three integers out.
					std::stringstream ss {line};
					int tmp_int {};
					while (ss >> tmp_int)
					{
						dims.push_back(tmp_int);	
					}

					// If we have anything other than 4 values something
					// went wrong.
					if (dims.size() != 4)
					{
						std::cerr << "Error! More than 4" 
							<< " dimensions were read in from " << filename 
							<< ". Check that this file is being generated " 
							<< "correctly.\n";
						std::cerr << "dims: ";
						for (auto i : dims) std::cerr << i << " ";
						std::cerr << '\n';
					}
				
					// nvalues is the total number of points in the file. The
					// data is flattened using C-style indexing, so as a 1D 
					// vector there are t*x*y*z values to read in. 
					ndata = std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<int>());
					ndata_read = true;
					std::cout << "Reading " << ndata << " data values\n";

					// Set the capacity of the vector now that we know what it will be.
					data_flattened.reserve(ndata);
				}
				else
				{
					data_flattened.push_back(std::stod(line));

					// If count ends up being greater than all the values
					// we need to read in, something went wrong.
					if (count >= ndata)
					{
						std::cerr << "Error! There appears to be too many "
							<< "values in " << filename << ". Check that "
							<< "this file is being generated correctly.\n";
						std::cerr << "count = " << count << "dims = ";
						for (std::size_t i {}; i < dims.size()-1; ++i)
						{
							std::cerr << dims[i] << ", ";
						}
						std::cerr << dims[dims.size()-1] << "\n";
					}

					++count;
				}
			}
		}
		fstream.close();

		// Move the data into a 4D vector and return it
		return {data_flattened, dims[0], dims[1], dims[2], dims[3]};
	}

	std::vector<std::string> load_interp_settings()
	{
		// Filename is treated as a constant.
		std::string filename {"bkg_from_pgkyl_interp_settings.csv"};
		std::ifstream fstream {filename};

		// The interpolation settings in a vector. FIrst element is
		// basis type, second is the poly order.
		std::vector<std::string> interp_settings {};

		// Go through one line at a time
		std::string line {};
		while(std::getline(fstream, line))
		{
			// If line starts with #, then it's a comment and ignore
			if (!line.starts_with("#"))
			{
				interp_settings.push_back(line);
			}
		}

		return interp_settings;
	}

	// Function to read in Gkeyll data using a python interface to postgkyl
	// via read_gkyl.py. This produces the following csv files:
	//   bkg_from_pgkyl_times.csv : The time for each frame
	//   bkg_from_pgkyl_grid.csv : Nodes for grid
	//   bkg_from_pgkyl_density.csv : Density arrays for all frames
	//   bkg_from_pgkyl_temperature.csv : Temperature arrays for all frames 
	// The data is loaded and placed into gkyl_data accordingly.
	void read_data_pgkyl(const std::string& species, 
		const std::string& data_type,
		Vectors::Vector4D& gkyl_data)
	{
		// Load into local variables so code is easier to read.
		const std::string gkyl_dir {Input::get_opt_str(Input::gkyl_dir)};
		const std::string gkyl_casename {Input::get_opt_str(Input::gkyl_casename)};
		const int gkyl_frame_start {Input::get_opt_int(Input::gkyl_frame_start)};
		const int gkyl_frame_end {Input::get_opt_int(Input::gkyl_frame_end)};
		const std::string gkyl_file_type {Input::get_opt_str(Input::gkyl_file_type)};

		// Ensure READ_GKYL_PY is defined, returning it if so.
		std::string read_gkyl_py {get_read_gkyl_py()};

		// Assemble command to call read_gkyl.py with stringstream	
		std::stringstream read_gkyl_cmd_ss {};
		read_gkyl_cmd_ss << "python " << read_gkyl_py 
			<< " --gkyl_dir=" << gkyl_dir
			<< " --gkyl_case_name=" << gkyl_casename
			<< " --gkyl_file_type=" << gkyl_file_type
			<< " --gkyl_frame_start=" << gkyl_frame_start
			<< " --gkyl_frame_end=" << gkyl_frame_end
			<< " --gkyl_species=" << species
			<< " --gkyl_data_type=" << data_type;

		// The bmag files don't include the needed DG interpolation data,
		// so we need to pass those in manually. To avoid user-error, we
		// will just rely on any of the other data file being created first,
		// where they will create a _interp_settings.csv file with the 
		// basis type and poly order that we need to pass in.
		if (data_type == "magnetic_field")
		{
			std::vector<std::string> interp_settings = load_interp_settings();
			read_gkyl_cmd_ss << " --gkyl_basis_type=" << interp_settings[0]
				<< " --gkyl_poly_order=" << interp_settings[1];
		}

		// Execute the command to save the files
		std::string tmp_str {read_gkyl_cmd_ss.str()};
		const char* read_gkyl_cmd {tmp_str.c_str()};
		std::cout << read_gkyl_cmd << '\n';
		[[maybe_unused]] int sys_result {system(read_gkyl_cmd)};

		// Load the times for each frame in a vector<double>
		std::vector<double> times {load_times()};

		// Load the x, y, z grids each as their own vector<double>
		auto [grid_x, grid_y, grid_z] = load_grid();

		// Load the values at each t, x, y, z into a Vector4D. We use
		// a tmp variable here instead of passing into move_into_data
		// directly because move_into_data uses move semantics, which
		// cannot be used on a const lvalue reference. So we create
		// the tmp rvalue and then pass it in, where the data is then
		// moved from tmp into gkyl_data. 
		std::cout << "Loading " << data_type << "...\n";
		Vectors::Vector4D tmp_data {load_values(data_type)};
		gkyl_data.move_into_data(tmp_data);
	}
/*
	// General function to read in data from Gkeyll into the relevant
	// vector specified by gkyl_data.
	void read_data(const std::string& species, const std::string& ftype,
		std::vector<Vectors::Vector3D>& gkyl_data, int comp)
	{
		// Load into local variables so code is easier to read.
		int gkyl_frame_start {Input::get_opt_int(Input::gkyl_frame_start)};
		int gkyl_frame_end {Input::get_opt_int(Input::gkyl_frame_end)};
		std::string gkyl_file_type {Input::get_opt_str(Input::gkyl_file_type)};
		
		// Loop through the frames.
		for (int f {gkyl_frame_start}; f <= gkyl_frame_end; ++f)
		{
			// Binary (.gkyl) type of file read in.
			if (gkyl_file_type == "binary")
			{
				// Assemble string of path to file.
				std::string path {assemble_path(species, ftype, f, ".gkyl")};
	
				// Read in a binary file.
				double time {};
				Vectors::Vector4D data {};
				std::tie(time, data) = GkylBinary::load_frame(path);
				std::cout << "Time = " << time << '\n';
				
				// The density is the first component in the last dimension.
				gkyl_data[f - gkyl_frame_start] = data.slice_dim4(comp);
			}
		}
	}
*/
	// Read electron density into gkyl_ne.
	void read_elec_density()
	{
		// Load into local variables so code is easier to read.
		std::string gkyl_elec_name {Input::get_opt_str(Input::gkyl_elec_name)};

		// Call pgkyl to load data into gkyl_ne
		read_data_pgkyl(gkyl_elec_name, "density", gkyl_ne);
	}

	// Read electron temperature into gkyl_te.
	void read_elec_temperature()
	{
		// Load into local variables so code is easier to read.
		std::string gkyl_elec_name {Input::get_opt_str(Input::gkyl_elec_name)};

		// Call pgkyl to load data into gkyl_te
		read_data_pgkyl(gkyl_elec_name, "temperature", gkyl_te);
	}

	// Read ion temperature into gkyl_ti.
	void read_ion_temperature()
	{
		// Load into local variables so code is easier to read.
		std::string gkyl_ion_name {Input::get_opt_str(Input::gkyl_ion_name)};

		// Call pgkyl to load data into gkyl_ti
		read_data_pgkyl(gkyl_ion_name, "temperature", gkyl_ti);
	}

	// Read plasma potential into gkyl_vp.
	void read_potential()
	{
		using namespace std::string_literals;

		// Call pgkyl to load data into gkyl_vp. There is no species
		// name with the potential file so just putting a placeholder
		// string. The script handles things.
		read_data_pgkyl("null"s, "potential", gkyl_vp);
	}

	// Read magnetic field into gkyl_b.
	void read_magnetic_field()
	{
		using namespace std::string_literals;

		// Call pgkyl to load data into gkyl_b. There is no species
		// name with the potential file so just putting a placeholder
		// string. The script handles things.
		read_data_pgkyl("null"s, "magnetic_field", gkyl_b);
	}

	Background::Background create_bkg()
	{
		Background::Background bkg {};

		// Move each of the vectors we've created into our bkg object
		bkg.move_into_ne(gkyl_ne);
		bkg.move_into_te(gkyl_te);
		bkg.move_into_ti(gkyl_ti);
		bkg.move_into_vp(gkyl_vp);
		bkg.move_into_b(gkyl_b);
		
		// Okay to return by value since C++11 uses move semantics here.
		return bkg;
	}
}
