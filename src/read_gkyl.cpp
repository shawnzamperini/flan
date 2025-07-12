/**
* @file read_gkyl.cpp
*
* @brief Routines handling reading in a plasma background from Gkeyll
*/

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>  // For std::multiplies
#include <iostream>
#include <numeric>     // For std::accumulate
#include <omp.h>
#include <string>
#include <sstream>
#include <type_traits> // For std::is_same
#include <tuple>

#include "background.h"
#include "constants.h"
#include "flan_types.h"
#include "options.h"
#include "read_gkyl.h"
#include "vectors.h"

namespace Gkyl
{
	// Entry point for reading Gkeyll data into Flan. 
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
		std::vector<BkgFPType> gkyl_times {};
		std::vector<BkgFPType> gkyl_x {};  // Cell centers
		std::vector<BkgFPType> gkyl_y {};
		std::vector<BkgFPType> gkyl_z {};
		std::vector<BkgFPType> gkyl_grid_x {};  // Grid edges
		std::vector<BkgFPType> gkyl_grid_y {};
		std::vector<BkgFPType> gkyl_grid_z {};
		Vectors::Vector3D<BkgFPType> gkyl_X {}; // Cell centers in physical space
		Vectors::Vector3D<BkgFPType> gkyl_Y {};
		Vectors::Vector3D<BkgFPType> gkyl_Z {};
		Vectors::Vector3D<BkgFPType> gkyl_grid_X {}; // Grid edges in physical 
		Vectors::Vector3D<BkgFPType> gkyl_grid_Y {}; // space
		Vectors::Vector3D<BkgFPType> gkyl_grid_Z {};
		Vectors::Vector4D<BkgFPType> gkyl_gij_00 {}; // Metric coefficients
		Vectors::Vector4D<BkgFPType> gkyl_gij_01 {}; 
		Vectors::Vector4D<BkgFPType> gkyl_gij_02 {};
		Vectors::Vector4D<BkgFPType> gkyl_gij_11 {};
		Vectors::Vector4D<BkgFPType> gkyl_gij_12 {};
		Vectors::Vector4D<BkgFPType> gkyl_gij_22 {};
		Vectors::Vector4D<BkgFPType> gkyl_J {};  // Jacobian
		Vectors::Vector4D<BkgFPType> gkyl_ne {};
		Vectors::Vector4D<BkgFPType> gkyl_te {};
		Vectors::Vector4D<BkgFPType> gkyl_ti {};
		Vectors::Vector4D<BkgFPType> gkyl_viperp_sq {};  // Perpendicular ion velocity
		Vectors::Vector4D<BkgFPType> gkyl_vp {};
		Vectors::Vector4D<BkgFPType> gkyl_bX {}; // Magnetic field components
		Vectors::Vector4D<BkgFPType> gkyl_bY {}; // in physical space
		Vectors::Vector4D<BkgFPType> gkyl_bZ {}; 
		Vectors::Vector4D<BkgFPType> gkyl_gradbX {}; // Magnetic field gradient
		Vectors::Vector4D<BkgFPType> gkyl_gradbY {}; // components in physical
		Vectors::Vector4D<BkgFPType> gkyl_gradbZ {}; // space
		Vectors::Vector4D<BkgFPType> gkyl_eX {}; // Electric field components
		Vectors::Vector4D<BkgFPType> gkyl_eY {}; // in physical space
		Vectors::Vector4D<BkgFPType> gkyl_eZ {};
		Vectors::Vector4D<BkgFPType> gkyl_uz {}; // Parallel (to z) ion flow
		Vectors::Vector4D<BkgFPType> gkyl_uX {}; // Flow components
		Vectors::Vector4D<BkgFPType> gkyl_uY {}; // in physical space
		Vectors::Vector4D<BkgFPType> gkyl_uZ {};

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
		//  grid_data[10] = gkyl_grid_X
		//  grid_data[11] = gkyl_grid_Y
		//  grid_data[12] = gkyl_grid_Z
		grid_data_t grid_data {gkyl_times, gkyl_x, gkyl_y, gkyl_z, 
			gkyl_grid_x, gkyl_grid_y, gkyl_grid_z, gkyl_X, gkyl_Y, gkyl_Z,
			gkyl_grid_X, gkyl_grid_Y, gkyl_grid_Z};

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
		read_magnetic_field(grid_data, gkyl_bX, gkyl_bY, gkyl_bZ, opts);

		// Calculate the Cartesian (X,Y,Z) coordinate of each cell center
		std::cout << "  - X,Y,Z coordinates\n";
		calc_cell_XYZ_centers(grid_data, opts);
		calc_cell_XYZ_edges(grid_data, opts);

		// Write out the (X,Y,Z) coordinates so that we can load them in python
		// to take advantage of scipy library in calculating gradients on
		// irregular grids (i.e., for the electric field). 
		write_XYZ(grid_data, opts);

		// Read in Jacobian
		std::cout << "  - Jacobian\n";
		read_jacobian(grid_data, gkyl_J, opts);

		// Read in metric coefficients. Each has its own vector.
		std::cout << "  - Metric coefficients\n";
		read_gij(grid_data, gkyl_gij_00, opts, "00");
		read_gij(grid_data, gkyl_gij_01, opts, "01");
		read_gij(grid_data, gkyl_gij_02, opts, "02");
		read_gij(grid_data, gkyl_gij_11, opts, "11");
		read_gij(grid_data, gkyl_gij_12, opts, "12");
		read_gij(grid_data, gkyl_gij_22, opts, "22");

		// At this point calculate gradients (E, gradB)
		std::cout << "Calculating gradient fields...\n";
		calc_gradients();

		// Read in electric field
		std::cout << "  - Electric field\n";
		read_elec_field(grid_data, gkyl_eX, gkyl_eY, gkyl_eZ);

		// These arrays are only needed if the gkyl moments are BiMaxwellian
		// (that is, it has vperp and vpar). We use vperp and gradB in the
		// background flow calculation, which are used in the collision
		// model. If you don't need collisions, then the other moment file
		// types are fine and this section is safely skipped.
		if (opts.gkyl_moment_type() == "bimaxwellian")
		{
			std::cout << "  - Ion parallel flow\n";
			read_par_flow(grid_data, gkyl_uz, opts);

			std::cout << "  - Ion perpendicular thermal velocity\n";
			read_ion_vperp_sq(grid_data, gkyl_viperp_sq, opts);

			// Read in magnetic field gradient. Note these are just needed to
			// calculate the gradB drift of the plasma flow, it is not needed
			// beyond that and thus doesn not go into bkg below.
			std::cout << "  - Magnetic field gradient\n";
			read_gradb(grid_data, gkyl_gradbX, gkyl_gradbY, gkyl_gradbZ);

			// Calculate the background flows in each Cartesian direction
			calc_bkg_flow(grid_data, gkyl_uz, gkyl_uX, 
				gkyl_uY, gkyl_uZ, gkyl_bX, gkyl_bY, gkyl_bZ, gkyl_eX, 
				gkyl_eY, gkyl_eZ, gkyl_gradbX, gkyl_gradbY, gkyl_gradbZ, 
				gkyl_viperp_sq, opts);
		}

		// With all our arrays assembled, encapsulate them into a Background
		// class object with move semantics then return.
		Background::Background bkg {};

		// Move each of the vectors we've created into our bkg object
		bkg.move_into_times(gkyl_times);
		bkg.move_into_x(gkyl_x);
		bkg.move_into_y(gkyl_y);
		bkg.move_into_z(gkyl_z);
		bkg.move_into_grid_x(gkyl_grid_x);
		bkg.move_into_grid_y(gkyl_grid_y);
		bkg.move_into_grid_z(gkyl_grid_z);
		bkg.move_into_ne(gkyl_ne);
		bkg.move_into_te(gkyl_te);
		bkg.move_into_ti(gkyl_ti);
		bkg.move_into_vp(gkyl_vp);
		bkg.move_into_bX(gkyl_bX);
		bkg.move_into_bY(gkyl_bY);
		bkg.move_into_bZ(gkyl_bZ);
		bkg.move_into_eX(gkyl_eX);
		bkg.move_into_eY(gkyl_eY);
		bkg.move_into_eZ(gkyl_eZ);
		bkg.move_into_uX(gkyl_uX);
		bkg.move_into_uY(gkyl_uY);
		bkg.move_into_uZ(gkyl_uZ);
		bkg.move_into_X(gkyl_X);
		bkg.move_into_Y(gkyl_Y);
		bkg.move_into_Z(gkyl_Z);
		bkg.move_into_grid_X(gkyl_grid_X);
		bkg.move_into_grid_Y(gkyl_grid_Y);
		bkg.move_into_grid_Z(gkyl_grid_Z);

		// Special case. gkyl_J is a 4D vector so we didn't have to write a
		// whole new function just for a 3D vector - less maintainable code. 
		// So we just need to index any random timeslice as a Vector3D, pass
		// it in. 
		Vectors::Vector3D gkyl_J_3D {gkyl_J.slice_dim1(0)};
		bkg.move_into_J(gkyl_J_3D);

		// Likewise for the 6 metric coefficient arrays
		Vectors::Vector3D gkyl_gij_00_3D {gkyl_gij_00.slice_dim1(0)};
		bkg.move_into_gij_00(gkyl_gij_00_3D);
		Vectors::Vector3D gkyl_gij_01_3D {gkyl_gij_01.slice_dim1(0)};
		bkg.move_into_gij_01(gkyl_gij_01_3D);
		Vectors::Vector3D gkyl_gij_02_3D {gkyl_gij_02.slice_dim1(0)};
		bkg.move_into_gij_02(gkyl_gij_02_3D);
		Vectors::Vector3D gkyl_gij_11_3D {gkyl_gij_11.slice_dim1(0)};
		bkg.move_into_gij_11(gkyl_gij_11_3D);
		Vectors::Vector3D gkyl_gij_12_3D {gkyl_gij_12.slice_dim1(0)};
		bkg.move_into_gij_12(gkyl_gij_12_3D);
		Vectors::Vector3D gkyl_gij_22_3D {gkyl_gij_22.slice_dim1(0)};
		bkg.move_into_gij_22(gkyl_gij_22_3D);

		// Okay to return by value since C++11 uses move semantics here.
		return bkg;
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
	std::vector<BkgFPType> load_times()
	{
		// Filename is treated as a constant.
		std::string filename {"bkg_from_pgkyl_times.csv"};
		std::ifstream fstream {filename};

		// Go through one line at a time
		std::string line {};
		std::vector<BkgFPType> times {};
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

	// Reads in the x, y, z grid nodes using pgkyl.
	std::tuple<std::vector<BkgFPType>, std::vector<BkgFPType>, std::vector<BkgFPType>> 
		load_grid()
	{
		// Filename is treated as a constant.
		std::string filename {"bkg_from_pgkyl_grid.csv"};
		std::ifstream fstream {filename};

		// Define variables needed for the following loop
		std::string line {};
		std::vector<BkgFPType> grid_x {};
		std::vector<BkgFPType> grid_y {};
		std::vector<BkgFPType> grid_z {};
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
					BkgFPType grid_value {static_cast<BkgFPType>(
						std::stod(line))};
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
	template <typename T>
	Vectors::Vector4D<T> load_values(const std::string& data_type)
	{
		// Filename is treated as a constant.
		std::string filename {"bkg_from_pgkyl_" + data_type + ".csv"};
		std::ifstream fstream {filename};

		// Define variables needed for the following loop
		std::string line {};
		std::vector<BkgFPType> data_flattened {};
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

	// Load file with the interpolation setting used by Gkeyll
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

	// Calculate cell centers for the given grid points.
	std::vector<BkgFPType> cell_centers(std::vector<BkgFPType>& grid)
	{

		std::vector<BkgFPType> centers (std::ssize(grid) - 1);
		for (std::size_t i {}; i < grid.size() - 1; ++i)
		{
			centers[i] = grid[i] + (grid[i+1] - grid[i]) / 2.0;
		}
		return centers;
	}

	// Function to read in Gkeyll data using a python interface to postgkyl
	// via read_gkyl.py. 
	template <typename T>
	void read_data_pgkyl(const std::string& species, 
		const std::string& data_type, grid_data_t& grid_data,
		Vectors::Vector4D<T>& gkyl_data, const Options::Options& opts, 
		const double species_mass_amu, const bool force_load)
	{
		// Check if file already exists and load it to save time. Can force
		// loading if we want to circumvent this. Account for the fact that
		// some files have the species name in them.
		std::string filename {};
		if (data_type == "density" || data_type == "temperature")
		{
			filename = "bkg_from_pgkyl_" + species + "_" + data_type + ".csv";
		}
		else
		{
			filename = "bkg_from_pgkyl_" + data_type + ".csv";
		}

		if (std::filesystem::exists(filename) && !force_load)
		{
			//std::cout << "Previous file located\n";
			if (data_type == "density" || data_type == "temperature")
			{
				Vectors::Vector4D<T> tmp_data {
					load_values<T>(species + "_" + data_type)};
				gkyl_data.move_into_data(tmp_data);
			}
			else
			{
				Vectors::Vector4D<T> tmp_data {load_values<T>(data_type)};
				gkyl_data.move_into_data(tmp_data);
			}
			return;
		}

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
			<< " --gkyl_data_type=" << data_type
			<< " --gkyl_moment_file_type=" << opts.gkyl_moment_type();

		// The bmag files don't include the needed DG interpolation data,
		// so we need to pass those in manually. To avoid user-error, we
		// will just rely on any of the other data file being created first,
		// where they will create a _interp_settings.csv file with the 
		// basis type and poly order that we need to pass in.
		if (data_type == "magnetic_unit_X" || 
			data_type == "magnetic_unit_Y" || 
			data_type == "magnetic_unit_Z" || 
			data_type == "magnetic_magnitude" || 
			data_type == "jacobian" ||
			data_type == "metric_coeff00" ||
			data_type == "metric_coeff01" ||
			data_type == "metric_coeff02" ||
			data_type == "metric_coeff11" ||
			data_type == "metric_coeff12" ||
			data_type == "metric_coeff22")
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
		if (data_type == "density" || data_type == "temperature")
		{
			// Some data is for specific species, which is going to be saved
			// as [species]_[data_type].
			Vectors::Vector4D<T> tmp_data {
				load_values<T>(species + "_" + data_type)};
			gkyl_data.move_into_data(tmp_data);
		}
		else
		{
			Vectors::Vector4D<T> tmp_data {load_values<T>(data_type)};
			gkyl_data.move_into_data(tmp_data);
		}
	}

	template <typename T>
	void read_elec_density(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_ne, const Options::Options& opts)
	{
		// Call pgkyl to load data into gkyl_ne. We force load the density
		// since it gets called first, and we need at least one load to 
		// load all the grid data.
		read_data_pgkyl(opts.gkyl_elec_name(), "density", grid_data, gkyl_ne, 
			opts, 0, true);
	}

	template <typename T>
	void read_elec_temperature(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_te, const Options::Options& opts)
	{
		// Call pgkyl to load data into gkyl_te. Need to pass in the species
		// mass so it can calculate temperature correctly.
		read_data_pgkyl(opts.gkyl_elec_name(), "temperature", grid_data,
			gkyl_te, opts, opts.gkyl_elec_mass_amu());
	}

	template <typename T>
	void read_ion_temperature(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_ti, const Options::Options& opts)
	{
		// Call pgkyl to load data into gkyl_ti. Need to pass in the species
		// mass so it can calculate temperature correctly.
		read_data_pgkyl(opts.gkyl_ion_name(), "temperature", grid_data, 
			gkyl_ti, opts, opts.gkyl_ion_mass_amu());
	}

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

	template <typename T>
	void read_magnetic_field(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_bX, Vectors::Vector4D<T>& gkyl_bY,
		Vectors::Vector4D<T>& gkyl_bZ, const Options::Options& opts)
	{
		using namespace std::string_literals;

		// There is no species name with the magnetic field files so just 
		// putting a placeholder string. The script handles things. 
		// These are actually the unit vector components of the magnetic field
		// so call once for each component. We multiply by the magnitude
		// later.
		read_data_pgkyl("null"s, "magnetic_unit_X", grid_data, gkyl_bX, opts);
		read_data_pgkyl("null"s, "magnetic_unit_Y", grid_data, gkyl_bY, opts);
		read_data_pgkyl("null"s, "magnetic_unit_Z", grid_data, gkyl_bZ, opts);

		// Load the magnetic field magnitude at each location
		Vectors::Vector4D<BkgFPType> gkyl_bmag {}; 
		read_data_pgkyl("null"s, "magnetic_magnitude", grid_data, gkyl_bmag, 
			opts);

		// The unit vector is actually a bit finnicky on the gkyl side of
		// things. It's not actually guaranteed to be 1.0 magnitude at the
		// cell centers which I guess is where it's calculated. So right now we
		// just hope the relative sizes of each component are mostly correct,
		// renormalize to 1.0 unit vector magnitude before multiplying by bmag.
		for (int i {}; i < std::ssize(std::get<0>(grid_data)); ++i)  // t
		{
		for (int j {}; j < std::ssize(std::get<1>(grid_data)); ++j)  // x
		{
		for (int k {}; k < std::ssize(std::get<2>(grid_data)); ++k)  // y
		{
		for (int l {}; l < std::ssize(std::get<3>(grid_data)); ++l)  // z
		{
			// Unit vector magnitude. Ideally 1.0, but life is not always so 
			// gentle. If the unit vector is just straight up wrong, then this
			// throws off the entire simulation. Need to asks a Gkeyll
			// developer if it's going to be correct or not.
			double unit_mag {gkyl_bX(i,j,k,l) * gkyl_bX(i,j,k,l)
				+ gkyl_bY(i,j,k,l) * gkyl_bY(i,j,k,l)
				+ gkyl_bZ(i,j,k,l) * gkyl_bZ(i,j,k,l)};
			
			// Renormalize each component, and multiply by the magnetic
			// field magnitude for each component.
			gkyl_bX(i,j,k,l) = gkyl_bX(i,j,k,l) * gkyl_bmag(i,j,k,l) 
				/ unit_mag;
			gkyl_bY(i,j,k,l) = gkyl_bY(i,j,k,l) * gkyl_bmag(i,j,k,l)
				/ unit_mag;
			gkyl_bZ(i,j,k,l) = gkyl_bZ(i,j,k,l) * gkyl_bmag(i,j,k,l) 
				/ unit_mag;
		}
		}
		}
		}
	}

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

	template <typename T>
	void read_gij(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_gij, const Options::Options& opts,
		const std::string& idx)
	{
		using namespace std::string_literals;

		// Assign the correct string for the index needed. idx is one of
		// 00, 01, 02, 11, 12, 22.
		std::string metric_coeff_str {"metric_coeff"};
		metric_coeff_str += idx;

		read_data_pgkyl("null"s, metric_coeff_str, grid_data, gkyl_gij, opts);
	}

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
	
	void calc_gradients()
	{
		// See if a file already exists and load it. If not, calculate it.
		// Let the user know one has been loaded, just in case they want Flan
		// to recalculate the values.
		std::string eX_fname {"bkg_from_pgkyl_elec_field_X.csv"};
		std::string eY_fname {"bkg_from_pgkyl_elec_field_Y.csv"};
		std::string eZ_fname {"bkg_from_pgkyl_elec_field_Z.csv"};
		std::string gradbX_fname {"bkg_from_pgkyl_gradb_X.csv"};
		std::string gradbY_fname {"bkg_from_pgkyl_gradb_Y.csv"};
		std::string gradbZ_fname {"bkg_from_pgkyl_gradb_Z.csv"};
		if ((std::filesystem::exists(eX_fname) && 
			std::filesystem::exists(eY_fname) &&
			std::filesystem::exists(eZ_fname) &&
			std::filesystem::exists(gradbX_fname) &&
			std::filesystem::exists(gradbY_fname) &&
			std::filesystem::exists(gradbZ_fname)))
		{
			//std::cout << "Previous file located\n";
		}
		else
		{
			// Execute python script to calculate electric field. First load
			// full path to python script.
			std::string calc_elec_field_py {get_calc_elec_field_py()};

			// Number of processes used is tied to OpenMP. We are not actually
			// using OpenMP here, we just need to see how many thread are 
			// available to us so we can pass it to the python script below.
			int num_processes {};
			#pragma omp parallel
			{
				if (omp_get_thread_num() == 0) num_processes 
					= omp_get_num_threads();
			}

			// The convert to a command
			std::stringstream calc_elec_field_cmd_ss {};
			calc_elec_field_cmd_ss << "python " << calc_elec_field_py
				<< " --num-processes " << num_processes; 
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

		// Resize the arrays since we know by now how big they need to be
		int dim1 {static_cast<int>(std::ssize(std::get<1>(grid_data)))};
		int dim2 {static_cast<int>(std::ssize(std::get<2>(grid_data)))};
		int dim3 {static_cast<int>(std::ssize(std::get<3>(grid_data)))};
		std::get<7>(grid_data).resize(dim1, dim2, dim3); // gkyl_X
		std::get<8>(grid_data).resize(dim1, dim2, dim3); // gkyl_Y
		std::get<9>(grid_data).resize(dim1, dim2, dim3); // gkyl_Z

		// Loop through every x, y, z value to calculate each X, Y, Z
		for (int i {}; i < dim1; ++i)
		{
			for (int j {}; j < dim2; ++j)
			{
				for (int k {}; k < dim3; ++k)
				{
					// Get the three Cartesian coordinates and store within
					// our 3D vector of 3-tuples
					auto [X, Y, Z] = opts.mapc2p()(std::get<1>(grid_data)[i], 
						std::get<2>(grid_data)[j], std::get<3>(grid_data)[k]);
					std::get<7>(grid_data)(i,j,k) = X;
					std::get<8>(grid_data)(i,j,k) = Y;
					std::get<9>(grid_data)(i,j,k) = Z;
				}
			}
		}
	}

	void calc_cell_XYZ_edges(grid_data_t& grid_data, 
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
		//  grid_data[10] = gkyl_grid_X
		//  grid_data[11] = gkyl_grid_Y
		//  grid_data[12] = gkyl_grid_Z

		// Resize the arrays since we know by now how big they need to be
		int dim1 {static_cast<int>(std::ssize(std::get<4>(grid_data)))};
		int dim2 {static_cast<int>(std::ssize(std::get<5>(grid_data)))};
		int dim3 {static_cast<int>(std::ssize(std::get<6>(grid_data)))};
		std::get<10>(grid_data).resize(dim1, dim2, dim3); // gkyl_grid_X
		std::get<11>(grid_data).resize(dim1, dim2, dim3); // gkyl_grid_Y
		std::get<12>(grid_data).resize(dim1, dim2, dim3); // gkyl_grid_Z

		// Loop through every x, y, z value to calculate each X, Y, Z
		for (int i {}; i < dim1; ++i)
		{
			//std::cout << i+1 << "/" << dim1 << '\n';
			for (int j {}; j < dim2; ++j)
			{
				for (int k {}; k < dim3; ++k)
				{

					// Get the three Cartesian coordinates and store within
					// our 3D vector of 3-tuples
					auto [X, Y, Z] = opts.mapc2p()(std::get<4>(grid_data)[i], 
						std::get<5>(grid_data)[j], std::get<6>(grid_data)[k]);
					std::get<10>(grid_data)(i,j,k) = X;
					std::get<11>(grid_data)(i,j,k) = Y;
					std::get<12>(grid_data)(i,j,k) = Z;
					//std::cout << "X, Y, Z, R: " << X << ", " << Y << ", "
					//	<< Z << ", " << std::sqrt(X*X+Y*Y) << '\n';
				}
			}
		}
	}

	void read_gradient_data(grid_data_t& grid_data, 
		Vectors::Vector4D<BkgFPType>& dataX, 
		Vectors::Vector4D<BkgFPType>& dataY, 
		Vectors::Vector4D<BkgFPType>& dataZ,
		Vectors::Vector4D<BkgFPType>& gkyl_dataX, 
		Vectors::Vector4D<BkgFPType>& gkyl_dataY, 
		Vectors::Vector4D<BkgFPType>& gkyl_dataZ)
	{

		// Apply fix to y boundaries where the gradient generally messes up
		gradient_ybound_fix(grid_data, dataX);
		gradient_ybound_fix(grid_data, dataY);
		gradient_ybound_fix(grid_data, dataZ);

		// Move into arrays
		gkyl_dataX.move_into_data(dataX);
		gkyl_dataY.move_into_data(dataY);
		gkyl_dataZ.move_into_data(dataZ);
	}

	void gradient_ybound_fix(grid_data_t& grid_data, 
		Vectors::Vector4D<BkgFPType>& ecomp)
	{
		// The following comment is for the electric field, but the same
		// logic applies to any other gradient-calculated data (e.g., the
		// magnetic field gradient).
		//
		// We apply a fix here. The y boundaries are treated as periodic (in
		// Gkeyll this need not always be true!). At the edge of the boundaries
		// the gradient calculation for the electric field can get messed up
		// because it does not have as many cells to calculate a gradient from.
		// So we mitigate this by assigning the edge cells a weighted average
		// from the "bounding" values in the y direction.
		// Consider if y has 6 cells:
		//
		//    |-----|-----|-----|
		// iy 4     5     0     1
		//             ^
		//          end of grid here, periodic
		//          so we restart at iy=0
		//
		// The problem cells are 0 and 5. We (arbitrarily) assign those as:
		//  iy=0   Ex = (1/3) * Ex_4 + (2/3) * Ex_1
		//  iy=5   Ex = (2/3) * Ex_4 + (1/3) * Ex_1
		// So we are taking the bounding two "reliable" E values at iy=1,4 and
		// reassigning the values at iy=0,5 with the weighted averages since
		// it is not clear what else we can do in this situation.
		int leny {static_cast<int>(std::ssize(std::get<2>(grid_data)))};
		for (int i {}; i < std::ssize(std::get<0>(grid_data)); i++)  // t
		for (int j {}; j < std::ssize(std::get<1>(grid_data)); j++)  // x
		for (int l {}; l < std::ssize(std::get<3>(grid_data)); l++)  // z
		{
		{
		{
			// Just naming these variables according to the example above.
			BkgFPType e1 {ecomp(i,j,1,l)};
			BkgFPType e4 {ecomp(i,j,leny-2,l)};
			ecomp(i,j,0,l) = (1.0/3.0) * e4 + (2.0/3.0) * e1;
			ecomp(i,j,leny-1,l) = (2.0/3.0) * e4 + (1.0/3.0) * e1;
		}
		}
		}
	}

	void read_elec_field(grid_data_t& grid_data, 
		Vectors::Vector4D<BkgFPType>& gkyl_eX, 
		Vectors::Vector4D<BkgFPType>& gkyl_eY, 
		Vectors::Vector4D<BkgFPType>& gkyl_eZ)
	{
		// Ex	
		Vectors::Vector4D<BkgFPType> tmp_dataX
			{load_values<BkgFPType>("elec_field_X")};

		// Ey
		Vectors::Vector4D<BkgFPType> tmp_dataY
			{load_values<BkgFPType>("elec_field_Y")};

		// Ez
		Vectors::Vector4D<BkgFPType> tmp_dataZ
			{load_values<BkgFPType>("elec_field_Z")};

		// Call general fuction to handle calculated-gradient data.
		read_gradient_data(grid_data, tmp_dataX, tmp_dataY, tmp_dataZ, 
			gkyl_eX, gkyl_eY, gkyl_eZ);

		// Apply fix to y boundaries where the gradient generally messes up
		/*
		elec_field_ybound_fix(grid_data, tmp_dataX);
		elec_field_ybound_fix(grid_data, tmp_dataY);
		elec_field_ybound_fix(grid_data, tmp_dataZ);

		// Move into arrays
		gkyl_eX.move_into_data(tmp_dataX);
		gkyl_eY.move_into_data(tmp_dataY);
		gkyl_eZ.move_into_data(tmp_dataZ);
		*/
	}

	void read_gradb(grid_data_t& grid_data, 
		Vectors::Vector4D<BkgFPType>& gkyl_gradbX, 
		Vectors::Vector4D<BkgFPType>& gkyl_gradbY, 
		Vectors::Vector4D<BkgFPType>& gkyl_gradbZ)
	{
		// gradB_X
		Vectors::Vector4D<BkgFPType> tmp_dataX
			{load_values<BkgFPType>("gradb_X")};

		// gradB_Y
		Vectors::Vector4D<BkgFPType> tmp_dataY
			{load_values<BkgFPType>("gradb_Y")};

		// gradB_Z
		Vectors::Vector4D<BkgFPType> tmp_dataZ
			{load_values<BkgFPType>("gradb_Z")};

		// Call general fuction to handle calculated-gradient data.
		read_gradient_data(grid_data, tmp_dataX, tmp_dataY, tmp_dataZ, 
			gkyl_gradbX, gkyl_gradbY, gkyl_gradbZ);
	}


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

		// Open up file for writing
		std::ofstream XYZ_file {"bkg_from_pgkyl_cell_XYZ.csv"};	
		if (!XYZ_file)
		{
			std::cerr << "Error creating bkg_from_pgkyl_cell_XYZ.csv!\n";
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

	
	template <typename T>
	void read_par_flow(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_uz, const Options::Options& opts)
	{
		// Call pgkyl to load data into gkyl_uz. This loads just the parallel
		// to z flow. The cross-field flows are drifts and added elsewhere.
		read_data_pgkyl(opts.gkyl_ion_name(), "par_flow", grid_data,
			gkyl_uz, opts);
	}

	template <typename T>
	void read_ion_vperp_sq(grid_data_t& grid_data, 
		Vectors::Vector4D<T>& gkyl_viperp_sq, const Options::Options& opts)
	{
		// Call pgkyl to load data into gkyl_viperp2. This requires
		// BiMaxwellian moments and thus will fail if they are not there.
		read_data_pgkyl(opts.gkyl_ion_name(), "vperp_sq", grid_data,
			gkyl_viperp_sq, opts);
	}

	// Ugly ass function header
	void calc_bkg_flow(grid_data_t& grid_data, 
		Vectors::Vector4D<BkgFPType>& gkyl_uz, 
		Vectors::Vector4D<BkgFPType>& gkyl_uX, 
		Vectors::Vector4D<BkgFPType>& gkyl_uY, 
		Vectors::Vector4D<BkgFPType>& gkyl_uZ,
		Vectors::Vector4D<BkgFPType>& gkyl_bX, 
		Vectors::Vector4D<BkgFPType>& gkyl_bY,
		Vectors::Vector4D<BkgFPType>& gkyl_bZ, 
		Vectors::Vector4D<BkgFPType>& gkyl_eX, 
		Vectors::Vector4D<BkgFPType>& gkyl_eY, 
		Vectors::Vector4D<BkgFPType>& gkyl_eZ, 
		Vectors::Vector4D<BkgFPType>& gkyl_gradbX, 
		Vectors::Vector4D<BkgFPType>& gkyl_gradbY,
		Vectors::Vector4D<BkgFPType>& gkyl_gradbZ, 
		Vectors::Vector4D<BkgFPType>& gkyl_viperp_sq, 
		const Options::Options& opts)
	{
		// Resize the uX, uY, uZ vectors
		gkyl_uX.resize(gkyl_uz.get_dim1(), gkyl_uz.get_dim2(), 
			gkyl_uz.get_dim3(), gkyl_uz.get_dim4());
		gkyl_uY.resize(gkyl_uz.get_dim1(), gkyl_uz.get_dim2(), 
			gkyl_uz.get_dim3(), gkyl_uz.get_dim4());
		gkyl_uZ.resize(gkyl_uz.get_dim1(), gkyl_uz.get_dim2(), 
			gkyl_uz.get_dim3(), gkyl_uz.get_dim4());

		// We are calculating uX = (parallel flow projection) + 
		// (ExB drift) + (gradB drift) = uz * bX + ExB_X + gradB_X
		for (int i {}; i < std::ssize(std::get<0>(grid_data)); ++i)  // t
		{
		for (int j {}; j < std::ssize(std::get<1>(grid_data)); ++j)  // x
		{
		for (int k {}; k < std::ssize(std::get<2>(grid_data)); ++k)  // y
		{
		for (int l {}; l < std::ssize(std::get<3>(grid_data)); ++l)  // z
		{
			// Some variables for cleaner code
			double bX {gkyl_bX(i,j,k,l)};
			double bY {gkyl_bY(i,j,k,l)};
			double bZ {gkyl_bZ(i,j,k,l)};
			double gradbX {gkyl_gradbX(i,j,k,l)};
			double gradbY {gkyl_gradbY(i,j,k,l)};
			double gradbZ {gkyl_gradbZ(i,j,k,l)};
			double eX {gkyl_eX(i,j,k,l)};
			double eY {gkyl_eY(i,j,k,l)};
			double eZ {gkyl_eZ(i,j,k,l)};
			double viperp_sq {gkyl_viperp_sq(i,j,k,l)};

			// Magnitude of magnetic field, used a couple times here
			double bmag_sq {bX * bX + bY * bY + bZ * bZ};
			double bmag {std::sqrt(bmag_sq)};

			// Projection of uz onto Cartesian magnetic field components
			gkyl_uX(i,j,k,l) = gkyl_uz(i,j,k,l) * bX / bmag;
			gkyl_uY(i,j,k,l) = gkyl_uz(i,j,k,l) * bY / bmag;
			gkyl_uZ(i,j,k,l) = gkyl_uz(i,j,k,l) * bZ / bmag;

			// ExB drift component
			gkyl_uX(i,j,k,l) = gkyl_uX(i,j,k,l) + (eY * bZ - eZ * bY) / bmag_sq;
			gkyl_uY(i,j,k,l) = gkyl_uY(i,j,k,l) + (eZ * bX - eX * bZ) / bmag_sq;
			gkyl_uZ(i,j,k,l) = gkyl_uZ(i,j,k,l) + (eX * bY - eY * bX) / bmag_sq;

			// gradB drift component. Group terms out front into coefficient.
			// Assume deuterium plasma (q = 1 C)
			double ion_charge {-Constants::charge_e};
			double coef {opts.gkyl_ion_mass_amu() * Constants::amu_to_kg 
				* viperp_sq / (2.0 * ion_charge * bmag_sq * bmag)};
			gkyl_uX(i,j,k,l) = gkyl_uX(i,j,k,l) + coef * (bY * gradbZ 
				- bZ * gradbY);
			gkyl_uY(i,j,k,l) = gkyl_uY(i,j,k,l) + coef * (bZ * gradbX 
				- bX * gradbZ);
			gkyl_uZ(i,j,k,l) = gkyl_uZ(i,j,k,l) + coef * (bX * gradbY 
				- bY * gradbX);

		}
		}
		}
		}
	}
}
