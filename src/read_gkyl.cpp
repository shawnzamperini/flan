/**
* @file read_gkyl.cpp
*
* @brief Routines handling reading in a plasma background from Gkeyll
*/

#include <filesystem>
#include <fstream>
#include <functional>  // For std::multiplies
#include <iostream>
#include <numeric>     // For std::accumulate
#include <string>
#include <sstream>
#include <type_traits> // For std::is_same
#include <tuple>

#include "background.h"
#include "constants.h"
#include "options.h"
#include "read_gkyl.h"
#include "vectors.h"

namespace Gkyl
{
	/**
	* @brief Entry point for reading Gkeyll data into Flan. 
	*
	* @param opts Options object containing the controlling setting of the
	* simulation.
	*
	* @return Returns a filled in Background object
	*/
	Background::Background read_gkyl(const Options::Options& opts)
	{	
		using namespace std::string_literals;

		// Vectors to hold all the various array loaded from a Gkeyll case.
		// Time-dependent data is 4D, and has dimensions (t, x, y, z). Time-
		// independent data is 3D and has dimension (x, y, z). In both cases,
		// the arrays are indexed according to the Gkeyll computational
		// coordinates indices. 
		// I feel like the logic could've shaken out better here, but in the 
		// end the data in these arrays are moved into a Background object 
		// (not copied) and then go out of scope, so end of the day it's fine. 
		std::vector<double> gkyl_times {};
		std::vector<double> gkyl_x {};  // Cell centers
		std::vector<double> gkyl_y {};
		std::vector<double> gkyl_z {};
		std::vector<int> gkyl_xidx {};  // Cell center idxs mapped to the
		std::vector<int> gkyl_yidx {};  // (flattened) X,Y,Z data
		std::vector<int> gkyl_zidx {};  
		std::vector<double> gkyl_grid_x {};  // Grid edges
		std::vector<double> gkyl_grid_y {};
		std::vector<double> gkyl_grid_z {};
		Vectors::Vector3D<double> gkyl_X {};
		Vectors::Vector3D<double> gkyl_Y {};
		Vectors::Vector3D<double> gkyl_Z {};
		Vectors::Vector4D<double> gkyl_J {};
		Vectors::Vector4D<double> gkyl_ne {};
		Vectors::Vector4D<double> gkyl_te {};
		Vectors::Vector4D<double> gkyl_ti {};
		Vectors::Vector4D<double> gkyl_vp {};
		Vectors::Vector4D<double> gkyl_b {}; 
		Vectors::Vector4D<double> gkyl_eX {};
		Vectors::Vector4D<double> gkyl_eY {};
		Vectors::Vector4D<double> gkyl_eZ {};

		// These are passed to every function call, so save some time/space by 
		// zipping them up into a tuple and passing together.
		//  grid_data[0]  = gkyl_times
		//  grid_data[1]  = gkyl_x
		//  grid_data[2]  = gkyl_y
		//  grid_data[3]  = gkyl_z
		//  grid_data[4]  = gkyl_grid_x
		//  grid_data[5]  = gkyl_grid_y
		//  grid_data[6]  = gkyl_grid_z
		//  grid_data[7]  = gkyl_X
		//  grid_data[8]  = gkyl_Y
		//  grid_data[9]  = gkyl_Z
		//  grid_data[10] = gkyl_xidx
		//  grid_data[11] = gkyl_yidx
		//  grid_data[12] = gkyl_zidx
		grid_data_t grid_data {gkyl_times, gkyl_x, gkyl_y, gkyl_z, 
			gkyl_grid_x, gkyl_grid_y, gkyl_grid_z, gkyl_X, gkyl_Y, gkyl_Z,
			gkyl_xidx, gkyl_yidx, gkyl_zidx};

		// Load each needed dataset from Gkeyll
		std::cout << "Loading Gkeyll data...\n";
		std::cout << "  - Electron density\n";
		read_elec_density(grid_data, gkyl_ne, opts);
		std::cout << "  - Electron temperature\n";
		read_elec_temperature(grid_data, gkyl_te, opts);
		std::cout << "  - Ion temperature\n";
		read_ion_temperature(grid_data, gkyl_ti, opts);
		std::cout << "  - Plasma potential\n";
		read_potential(grid_data, gkyl_vp, opts);
		std::cout << "  - Magnetic field\n";
		read_magnetic_field(grid_data, gkyl_b, opts);

		// Calculate the Cartesian (X,Y,Z) coordinate of each cell center
		std::cout << "  - X,Y,Z coordinates\n";
		calc_cell_XYZ_centers(grid_data, opts);

		// Write out the (X,Y,Z) coordinates so that we can load them in python
		// to take advantage of scipy library in calculating gradients on
		// irregular grids (i.e., for the electric field). 
		write_XYZ(grid_data, opts);

		// Calculate the 3 electric field components from the potential
		calc_elec_field();

		// Then read each electric field component
		std::cout << "  - Electric field\n";
		read_elec_field(gkyl_eX, gkyl_eY, gkyl_eZ);

		// Read in Jacobian
		std::cout << "  - Jacobian\n";
		read_jacobian(grid_data, gkyl_J, opts);

		// With all our arrays assembled, encapsulate them into a Background
		// class object with move semantics then return.
		//Background::Background bkg {create_bkg()};
		Background::Background bkg {};

		// Move each of the vectors we've created into our bkg object
		bkg.move_into_times(gkyl_times);
		bkg.move_into_x(gkyl_x);
		bkg.move_into_y(gkyl_y);
		bkg.move_into_z(gkyl_z);
		bkg.move_into_xidx(gkyl_xidx);
		bkg.move_into_yidx(gkyl_yidx);
		bkg.move_into_zidx(gkyl_zidx);
		bkg.move_into_grid_x(gkyl_grid_x);
		bkg.move_into_grid_y(gkyl_grid_y);
		bkg.move_into_grid_z(gkyl_grid_z);
		bkg.move_into_ne(gkyl_ne);
		bkg.move_into_te(gkyl_te);
		bkg.move_into_ti(gkyl_ti);
		bkg.move_into_vp(gkyl_vp);
		bkg.move_into_b(gkyl_b);
		bkg.move_into_eX(gkyl_eX);
		bkg.move_into_eY(gkyl_eY);
		bkg.move_into_eZ(gkyl_eZ);
		bkg.move_into_X(gkyl_X);
		bkg.move_into_Y(gkyl_Y);
		bkg.move_into_Z(gkyl_Z);

		// Special case. gkyl_J is a 4D vector so we didn't have to write a
		// whole new function just for a 3D vector - less maintainable code. 
		// So we just need to index any random timeslice as a Vector3D, pass
		// it in. 
		Vectors::Vector3D gkyl_J_3D {gkyl_J.slice_dim1(0)};
		bkg.move_into_J(gkyl_J_3D);

		// Okay to return by value since C++11 uses move semantics here.
		return bkg;
	}
	
	/**
	* @brief Helper function to get Gkeyll file path
	*
	* Simple helper function to return a string of the full path to a Gkeyll
	* file using it's file name convention.
	*
	* @param species Name of the species used in the Gkeyll run
	* @param ftype The data type to load, e.g., M0 or prim_moms.
	* @param frame The frame number
	* @param extension The extension type, e.g.,\ .gkyl.
	*
	* @return Returns a string of the full file path.
	*/
	/*
	std::string assemble_path(const std::string& species, 
		const std::string& ftype, int frame, const std::string& extension,
		const Options::Options& opts) 
	{
		return opts.gkyl_dir() + "/" + opts.gkyl_casename() + "-" + species 
			+ "_" + ftype + "_" + std::to_string(frame) + extension;
	}
	*/

	/**
	* @brief Get path to python interface script to postgkyl, read_gkyl\.py.
	*
	* Function to return the full path to read_gkyl.py, which is used to
	* interface with postgkyl and create files that are easily read in by
	* Flan. This path is stored as an environment variable that is set
	* as part of the flan conda environment.
	*
	* @return Returns string containing the full path to read_gkyl.py.
	*/
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

	/**
	* @brief Reads in a 1D vector of times corresponding to each frame.
	*
	* @return Return a vector of the times.
	*/
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
					//std::cout << "Reading " << ntimes << " times\n";

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

	/**
	* @brief Reads in the x, y, z grid nodes using pgkyl.
	*
	* @return Returns a tuple of 3 vectors containing the grid edges - (x,y,z)
	*/
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
		[[maybe_unused]] int ngrid {};
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
					//std::cout << "Reading " << ngrid << " grid values\n";

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

	/**
	* @brief Read in data values using pgkyl, returning as a Vector4D.
	*
	* @param data_type String identifying which data file to load. This is one
	* of density, temperature, magnetic_field, potential or times.
	*
	* @return Returns a Vector4D object of the data (t,x,y,z).
	*/
	template <typename T>
	Vectors::Vector4D<T> load_values(const std::string& data_type)
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
					ndata = std::accumulate(dims.begin(), dims.end(), 1, 
						std::multiplies<int>());
					ndata_read = true;
					//std::cout << "Reading " << ndata << " data values\n";

					// Set the capacity of the vector now that we know what it 
					// will be.
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

	/**
	* @brief Load file with the interpolation setting used by Gkeyll
	*
	* This file is made automatically so it is available to other files which
	* do not have this data included with them (like the magnetic_field data).
	*
	* @returns Returns a vector of two elements: the basis type and the order.
	*/
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

	/**
	* @brief Calculate cell centers for the given grid points.
	*
	* Calculate the cell centers for the given grid points. This returns a
	* vector one shorter than the grid point. This only works for uniform
	* grids, just in case that isn't obvious throughout the whole code.
	*
	* @param grid The grid edges for a given dimension.
	* @return Returns a vector of the cell centers that is one less in length
	* than grid.
	*/
	std::vector<double> cell_centers(std::vector<double>& grid)
	{

		std::vector<double> centers (std::ssize(grid) - 1);
		for (std::size_t i {}; i < grid.size() - 1; ++i)
		{
			centers[i] = grid[i] + (grid[i+1] - grid[i]) / 2.0;
		}
		return centers;
	}

	/**
	* @brief Function to read in Gkeyll data using a python interface to postgkyl
	* via read_gkyl.py. 
	*
	* This produces the following csv files:
	*   bkg_from_pgkyl_times.csv : The time for each frame
	*   bkg_from_pgkyl_grid.csv : Nodes for grid
	*   bkg_from_pgkyl_density.csv : Density arrays for all frames
	*   bkg_from_pgkyl_temperature.csv : Temperature arrays for all frames 
	*   bkg_from_pgkyl_potential.csv : Plasma potential arrays for all frames 
	*   bkg_from_pgkyl_magnetic_field.csv : Magnetic field arrays for all frames 
	* The data is loaded and placed into gkyl_data accordingly.
	*
	* @param species Name of the species for the Gkeyll run
	* @param data_type Name of the data to load.\ Can be one of temperature,
	* density, potential or magnetic field. 
	* @param gkyl_data One of the global (to the Gkyl namespace) vectors that
	* corredponds to the data_type being loaded.
	* @param opts Options object that contains all the controlling options
	* of the simulation.
	* @param species_mass_amu The mass of the species in amu.
	*/
	template <typename T>
	void read_data_pgkyl(const std::string& species, 
		const std::string& data_type, grid_data_t& grid_data,
		Vectors::Vector4D<T>& gkyl_data, const Options::Options& opts, 
		const double species_mass_amu)
	{
		// Ensure READ_GKYL_PY is defined, returning it if so.
		std::string read_gkyl_py {get_read_gkyl_py()};

		// Assemble command to call read_gkyl.py with stringstream	
		std::stringstream read_gkyl_cmd_ss {};
		read_gkyl_cmd_ss << "python " << read_gkyl_py 
			<< " --gkyl_dir=" << opts.gkyl_dir()
			<< " --gkyl_case_name=" << opts.gkyl_casename()
			<< " --gkyl_file_type=" << opts.gkyl_file_type()
			<< " --gkyl_frame_start=" << opts.gkyl_frame_start()
			<< " --gkyl_frame_end=" << opts.gkyl_frame_end()
			<< " --gkyl_species=" << species
			<< " --gkyl_species_mass_amu=" << species_mass_amu
			<< " --gkyl_data_type=" << data_type;

		// The bmag files don't include the needed DG interpolation data,
		// so we need to pass those in manually. To avoid user-error, we
		// will just rely on any of the other data file being created first,
		// where they will create a _interp_settings.csv file with the 
		// basis type and poly order that we need to pass in.
		if ((data_type == "magnetic_field") || (data_type == "jacobian"))
		{
			std::vector<std::string> interp_settings = load_interp_settings();
			read_gkyl_cmd_ss << " --gkyl_basis_type=" << interp_settings[0]
				<< " --gkyl_poly_order=" << interp_settings[1];
		}

		// Execute the command to save the files
		std::string tmp_str {read_gkyl_cmd_ss.str()};
		const char* read_gkyl_cmd {tmp_str.c_str()};
		//std::cout << read_gkyl_cmd << '\n';
		[[maybe_unused]] int sys_result {system(read_gkyl_cmd)};

		// Load the times for each frame in a vector<double>
		//gkyl_times = {load_times()};
		std::get<0>(grid_data) = {load_times()};  // gkyl_times

		// Load the x, y, z grids each as their own vector<double>
		//std::tie(gkyl_grid_x, gkyl_grid_y, gkyl_grid_z) = load_grid();
		std::tie(std::get<4>(grid_data), std::get<5>(grid_data), 
			std::get<6>(grid_data)) = load_grid();

		// Calculate the cell centers.
		//gkyl_x = cell_centers(gkyl_grid_x);
		//gkyl_y = cell_centers(gkyl_grid_y);
		//gkyl_z = cell_centers(gkyl_grid_z);
		std::get<1>(grid_data) = cell_centers(std::get<4>(grid_data));
		std::get<2>(grid_data) = cell_centers(std::get<5>(grid_data));
		std::get<3>(grid_data) = cell_centers(std::get<6>(grid_data));

		// Load the values at each t, x, y, z into a Vector4D. We use
		// a tmp variable here instead of passing into move_into_data
		// directly because move_into_data uses move semantics, which
		// cannot be used on a const lvalue reference. So we create
		// the tmp rvalue and then pass it in, where the data is then
		// moved from tmp into gkyl_data. 
		//std::cout << "Loading " << data_type << "...\n";
		Vectors::Vector4D<T> tmp_data {load_values<T>(data_type)};
		gkyl_data.move_into_data(tmp_data);
	}

	/**
	* @brief Read electron density into gkyl_ne.
	*/
	template <typename T>
	void read_elec_density(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_ne, const Options::Options& opts)
	{
		// Call pgkyl to load data into gkyl_ne
		read_data_pgkyl(opts.gkyl_elec_name(), "density", grid_data, gkyl_ne, 
			opts);
	}

	/**
	* @brief Read electron temperature into gkyl_te.
	*/
	template <typename T>
	void read_elec_temperature(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_te, const Options::Options& opts)
	{
		// Call pgkyl to load data into gkyl_te. Need to pass in the species
		// mass so it can calculate temperature correctly.
		read_data_pgkyl(opts.gkyl_elec_name(), "temperature", grid_data,
			gkyl_te, opts, opts.gkyl_elec_mass_amu());
	}

	/**
	* @brief Read ion temperature into gkyl_ti.
	*/
	template <typename T>
	void read_ion_temperature(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_ti, const Options::Options& opts)
	{
		// Call pgkyl to load data into gkyl_ti. Need to pass in the species
		// mass so it can calculate temperature correctly.
		read_data_pgkyl(opts.gkyl_ion_name(), "temperature", grid_data, 
			gkyl_ti, opts, opts.gkyl_ion_mass_amu());
	}

	/**
	* @brief Read plasma potential into gkyl_vp.
	*/
	template <typename T>
	void read_potential(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_vp, const Options::Options& opts)
	{
		using namespace std::string_literals;

		// Call pgkyl to load data into gkyl_vp. There is no species
		// name with the potential file so just putting a placeholder
		// string. The script handles things.
		read_data_pgkyl("null"s, "potential", grid_data, gkyl_vp, opts);
	}

	/**
	* @brief Read magnetic field into gkyl_b.
	*/
	template <typename T>
	void read_magnetic_field(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_b, const Options::Options& opts)
	{
		using namespace std::string_literals;

		// Call pgkyl to load data into gkyl_b. There is no species
		// name with the potential file so just putting a placeholder
		// string. The script handles things.
		read_data_pgkyl("null"s, "magnetic_field", grid_data, gkyl_b, opts);
	}

	/**
	* @brief Read Jacobian into gkyl_J.
	*/
	template <typename T>
	void read_jacobian(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_J, const Options::Options& opts)
	{
		using namespace std::string_literals;

		// Call pgkyl to load data into gkyl_b. There is no species
		// name with the potential file so just putting a placeholder
		// string. The script handles things.
		read_data_pgkyl("null"s, "jacobian", grid_data, gkyl_J, opts);
	}

	/**
	* @brief Get path to python script for calculating the electric field, 
	* calc_elec_field.py.
	*
	* Function to return the full path to calc_elec_field.py, which uses
	* the files cell_center_XYZ.csv and bkg_from_pgkyl_potential.csv to
	* calculate the electric field components from the potential gradient.
	* This path is stored as an environment variable that is set
	* as part of the flan conda environment.
	*
	* @return Returns string containing the full path to calc_elec_field.py.
	*/
	std::string get_calc_elec_field_py()
	{
		if (const char* calc_elec_field_py {std::getenv("CALC_ELEC_FIELD_PY")})
		{
			return static_cast<std::string>(calc_elec_field_py);
		}
		else
		{
			std::cerr << "Error! CALC_ELEC_FIELD_PY is not set. Please set "
				<< "this to the full path to the calc_elec_field.py script. "
				<< "E.g.:\n"
				<< "  CALC_ELEC_FIELD_PY=/path/to/flan/python/"
				<< "calc_elec_field.py\n";
			return std::string {};
		}
	}
	void calc_elec_field()
	{
		// See if a file already exists and load it. If not, calculate it.
		// Let the user know one has been loaded, just in case they want Flan
		// to recalculate the values.
		std::string eX_fname {"bkg_from_pgkyl_elec_field_X.csv"};
		std::string eY_fname {"bkg_from_pgkyl_elec_field_Y.csv"};
		std::string eZ_fname {"bkg_from_pgkyl_elec_field_Z.csv"};
		if ((std::filesystem::exists(eX_fname) && 
			std::filesystem::exists(eY_fname) &&
			std::filesystem::exists(eZ_fname)))
		{
			std::cout << "Electric field files located. Skipping "
				<< "calculation.\n";
		}
		else
		{
			// Execute python script to calculate electric field. First load
			// full path to python script.
			std::string calc_elec_field_py {get_calc_elec_field_py()};

			// The convert to a command
			std::stringstream calc_elec_field_cmd_ss {};
			calc_elec_field_cmd_ss << "python " << calc_elec_field_py; 
			std::string tmp_str {calc_elec_field_cmd_ss.str()};
			const char* calc_elec_field_cmd {tmp_str.c_str()};
				
			// And execute it
			[[maybe_unused]] int sys_result {system(calc_elec_field_cmd)};
		}
	}
	
	void calc_cell_XYZ_centers(grid_data_t& grid_data, 
		const Options::Options& opts)
	{
		//  grid_data[0]  = gkyl_times
		//  grid_data[1]  = gkyl_x
		//  grid_data[2]  = gkyl_y
		//  grid_data[3]  = gkyl_z
		//  grid_data[4]  = gkyl_grid_x
		//  grid_data[5]  = gkyl_grid_y
		//  grid_data[6]  = gkyl_grid_z
		//  grid_data[7]  = gkyl_X
		//  grid_data[8]  = gkyl_Y
		//  grid_data[9]  = gkyl_Z
		//  grid_data[10] = gkyl_xidx
		//  grid_data[11] = gkyl_yidx
		//  grid_data[12] = gkyl_zidx

		// Resize the arrays since we know by now how big they need to be
		//int dim1 {static_cast<int>(std::ssize(gkyl_x))};
		//int dim2 {static_cast<int>(std::ssize(gkyl_y))};
		//int dim3 {static_cast<int>(std::ssize(gkyl_z))};
		int dim1 {static_cast<int>(std::ssize(std::get<1>(grid_data)))};
		int dim2 {static_cast<int>(std::ssize(std::get<2>(grid_data)))};
		int dim3 {static_cast<int>(std::ssize(std::get<3>(grid_data)))};
		//gkyl_X.resize(dim1, dim2, dim3);
		//gkyl_Y.resize(dim1, dim2, dim3);
		//gkyl_Z.resize(dim1, dim2, dim3);
		std::get<7>(grid_data).resize(dim1, dim2, dim3); // gkyl_X
		std::get<8>(grid_data).resize(dim1, dim2, dim3); // gkyl_Y
		std::get<9>(grid_data).resize(dim1, dim2, dim3); // gkyl_Z
		//gkyl_xidx.resize(dim1 * dim2 * dim3);
		//gkyl_yidx.resize(dim1 * dim2 * dim3);
		//gkyl_zidx.resize(dim1 * dim2 * dim3);
		std::get<10>(grid_data).resize(dim1 * dim2 * dim3); // gkyl_xidx
		std::get<11>(grid_data).resize(dim1 * dim2 * dim3); // gkyl_yidx
		std::get<12>(grid_data).resize(dim1 * dim2 * dim3); // gkyl_zidx

		// Loop through every x, y, z value to calculate each X, Y, Z
		for (int i {}; i < std::ssize(std::get<1>(grid_data)); ++i)
		{
			for (int j {}; j < std::ssize(std::get<2>(grid_data)); ++j)
			{
				for (int k {}; k < std::ssize(std::get<3>(grid_data)); ++k)
				{
					// Get the three Cartesian coordinates and store within
					// our 3D vector of 3-tuples
					auto [X, Y, Z] = opts.mapc2p()(std::get<1>(grid_data)[i], 
						std::get<2>(grid_data)[j], std::get<3>(grid_data)[k]);
					std::get<7>(grid_data)(i,j,k) = X;
					std::get<8>(grid_data)(i,j,k) = Y;
					std::get<9>(grid_data)(i,j,k) = Z;

					// Store the indices in computational coordinates that this
					// physical coordinate cooresponds to. Can just use
					// calc_index from the 3D vector (each X, Y, Z calc_index
					// gives the same thing).
					std::get<10>(grid_data)[std::get<7>(grid_data)
						.calc_index(i,j,k)] = i;
					std::get<11>(grid_data)[std::get<8>(grid_data)
						.calc_index(i,j,k)] = j;
					std::get<12>(grid_data)[std::get<9>(grid_data)
						.calc_index(i,j,k)] = k;
				}
			}
		}
	}

	void read_elec_field(Vectors::Vector4D<double>& gkyl_eX, 
		Vectors::Vector4D<double>& gkyl_eY, Vectors::Vector4D<double>& gkyl_eZ)
	{
		// Ex	
		Vectors::Vector4D<double> tmp_dataX
			{load_values<double>("elec_field_X")};
		gkyl_eX.move_into_data(tmp_dataX);

		// Ey
		Vectors::Vector4D<double> tmp_dataY
			{load_values<double>("elec_field_Y")};
		gkyl_eY.move_into_data(tmp_dataY);

		// Ez
		Vectors::Vector4D<double> tmp_dataZ
			{load_values<double>("elec_field_Z")};
		gkyl_eZ.move_into_data(tmp_dataZ);
	}

	/**
	* @brief Write out Cartesian (X,Y,Z) coordinates of cell centers
	*
	* This function writes to a file of X,Y,Z coordinates,
	* so that we can read them into a python script and take advantage of
	* the scipy library to easily calculate a gradient on an irregular grid.
	*/
	void write_XYZ(grid_data_t& grid_data, const Options::Options& opts)
	{
		//  grid_data[0]  = gkyl_times
		//  grid_data[1]  = gkyl_x
		//  grid_data[2]  = gkyl_y
		//  grid_data[3]  = gkyl_z
		//  grid_data[4]  = gkyl_grid_x
		//  grid_data[5]  = gkyl_grid_y
		//  grid_data[6]  = gkyl_grid_z
		//  grid_data[7]  = gkyl_X
		//  grid_data[8]  = gkyl_Y
		//  grid_data[9]  = gkyl_Z
		//  grid_data[10] = gkyl_xidx
		//  grid_data[11] = gkyl_yidx
		//  grid_data[12] = gkyl_zidx

		// Open up file for writing
		std::ofstream XYZ_file {"cell_center_XYZ.csv"};	
		if (!XYZ_file)
		{
			std::cerr << "Error creating cell_center_XYZ.csv!\n";
		}

		// Print header info
		XYZ_file << "# This file contains the Cartesian X,Y,Z coordinates of\n"
			     << "# the Gkeyll grid at each cell center. They are\n"
				 << "# calculated using the mapc2p function passed in from\n"
				 << "# the input file. The first line in this file are the\n"
				 << "# dimensions of the data. The coordinates\n"
				 << "# are then printed one at a time with X Y Z on each \n"
				 << "# line.\n";

		// Print out shape of data
		XYZ_file << std::get<7>(grid_data).get_dim1() << " " 
			<< std::get<7>(grid_data).get_dim2() << " " 
			<< std::get<7>(grid_data).get_dim3() << '\n';

		// Print out the coordinates
		for (int i {}; i < std::ssize(std::get<1>(grid_data)); ++i)
		{
			for (int j {}; j < std::ssize(std::get<2>(grid_data)); ++j)
			{
				for (int k {}; k < std::ssize(std::get<3>(grid_data)); ++k)
				{
					XYZ_file << std::get<7>(grid_data)(i,j,k) << " " 
						<< std::get<8>(grid_data)(i,j,k) << " " 
						<< std::get<9>(grid_data)(i,j,k) << '\n';
				}
			}
		}
	}

	/**
	* @brief Move data from the global (to the Gkyl namespace) into a
	* Background object and return.
	*
	* @return Returns a filled-in Background object to perform an impurity
	* transport simulation on.
	*/
	/*
	Background::Background create_bkg()
	{
		Background::Background bkg {};

		// Move each of the vectors we've created into our bkg object
		bkg.move_into_times(gkyl_times);
		bkg.move_into_x(gkyl_x);
		bkg.move_into_y(gkyl_y);
		bkg.move_into_z(gkyl_z);
		bkg.move_into_xidx(gkyl_xidx);
		bkg.move_into_yidx(gkyl_yidx);
		bkg.move_into_zidx(gkyl_zidx);
		bkg.move_into_grid_x(gkyl_grid_x);
		bkg.move_into_grid_y(gkyl_grid_y);
		bkg.move_into_grid_z(gkyl_grid_z);
		bkg.move_into_ne(gkyl_ne);
		bkg.move_into_te(gkyl_te);
		bkg.move_into_ti(gkyl_ti);
		bkg.move_into_vp(gkyl_vp);
		bkg.move_into_b(gkyl_b);
		bkg.move_into_eX(gkyl_eX);
		bkg.move_into_eY(gkyl_eY);
		bkg.move_into_eZ(gkyl_eZ);
		bkg.move_into_X(gkyl_X);
		bkg.move_into_Y(gkyl_Y);
		bkg.move_into_Z(gkyl_Z);

		// Okay to return by value since C++11 uses move semantics here.
		return bkg;
	}
	*/
}
