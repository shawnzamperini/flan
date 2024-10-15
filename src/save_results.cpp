/**
* @file save_results.cpp
* @brief Routines for saving results of Flan simulation
*
* Right now only outputting everything into a single netCDF file is supported.
*/

#include <type_traits>
#include <iostream>
#include <string>
#include <vector>
#include <variant>

#include "background.h"
#include "save_results.h"
#include "netcdf"

namespace SaveResults
{
	void save_results(const std::string_view case_name, 
		const Background::Background& bkg, Impurity::Statistics& imp_stats)
	{
		// Right now this is just a redirect to the netcdf option. I built 
		// it this way to allow other save options in the future if desired.
		save_netcdf(case_name, bkg, imp_stats);
	
	}

	// Save a 1D vector into nc_file. dim must be matched by the programmer.
	template <typename T>
	void save_vector_1d(const netCDF::NcFile& nc_file, std::vector<T> vec, 
		const std::string& var_name, const netCDF::NcDim& dim, 
		const std::string& description, const std::string& units)
	{
		// Vector cannot be zero-length.
		if (vec.size() == 0)
		{
			std::cerr << "NetCDF Error! " << var_name << " is zero-length."
				<< " Skipping.\n";
		}
		else
		{
			// Create variable with netCDF type based on which template is being 
			// used.
			netCDF::NcVar var {};
			if constexpr (std::is_same_v<T, double>)
			{
				var = nc_file.addVar(var_name, netCDF::ncDouble, dim);
			}
			else if constexpr (std::is_same_v<T, int>)
			{
				var = nc_file.addVar(var_name, netCDF::ncInt, dim);
			}

			// Put vector data into it.
			var.putVar(vec.data());

			// Add description
			var.putAtt("description", description);

			// Add units
			var.putAtt("units", units);
		}
	}

	// Save a 4D vector into nc_file. dims must be matched by the programmer.
	template <typename T>
	void save_vector_4d(const netCDF::NcFile& nc_file, 
		const Vectors::Vector4D<T>& vec, 
		const std::string& var_name, const netCDF::NcDim& dim1,
		const netCDF::NcDim& dim2, const netCDF::NcDim& dim3, 
		const netCDF::NcDim& dim4, const std::string& description, 
		const std::string& units)
	{
		if (vec.get_data().size() == 0)
		{
			std::cerr << "NetCDF Error! " << var_name << " is zero-length."
				<< " Skipping.\n";
		}
		else
		{
			// Create a 4D variable	
			netCDF::NcVar var {nc_file.addVar(var_name, netCDF::ncDouble, 
					{dim1, dim2, dim3, dim4})};

			// We actually save it as a flattened array no matter the
			// dimension, which is how we represent the data internally 
			// already. So just throw that data in there.
			var.putVar(vec.get_data().data());

			// Add description
			var.putAtt("description", description);

			// Add units
			var.putAtt("units", units);
		}
	}

	void save_netcdf(const std::string_view case_name, 
		const Background::Background& bkg, Impurity::Statistics& imp_stats)
	{
		std::cout << "Saving netCDF file...\n";

		// Reusable string to pass in descriptions
		std::string desc {};

		// Open new netCDF file named after this case. replace tells netCDF
		// to overwrite the file if it exists.
		// NcFile's scope is private to whatever block it is defined in,
		// so we cannot pass it around functions (at least not like a normal
		// variable). The implementation details are such that it automatically
		// closes once it goes out of scope.
		std::string nc_filename {case_name};
		nc_filename = nc_filename + ".nc";
		netCDF::NcFile nc_file {nc_filename, netCDF::NcFile::replace};
		
		// Case name
		//nc_file.putAtt("case_name", case_name);

		// Details about the simulation (system, time run, time spent, etc.)
		// To-do

		// Save all input options used for this case
		// To-do

		// Dimensions for the arrays containing data at the cell centers. 
		// First dimension is time, then x, y, z.
		std::cout << "Creating dimensions...";
		netCDF::NcDim dim1 {nc_file.addDim("time", bkg.get_dim1())};
		netCDF::NcDim dim2 {nc_file.addDim("x", bkg.get_dim2())};
		netCDF::NcDim dim3 {nc_file.addDim("y", bkg.get_dim3())};
		netCDF::NcDim dim4 {nc_file.addDim("z", bkg.get_dim4())};
		std::cout << " done\n";

		// Dimensions for the grid points. They're just bigger than the
		// cell center lengths.
		netCDF::NcDim grid_dim2 {nc_file.addDim("grid x", bkg.get_dim2() + 1)};
		netCDF::NcDim grid_dim3 {nc_file.addDim("grid y", bkg.get_dim3() + 1)};
		netCDF::NcDim grid_dim4 {nc_file.addDim("grid z", bkg.get_dim4() + 1)};

		// Time
		desc = "time for each frame";
		save_vector_1d(nc_file, bkg.get_times(), "time", dim1, desc, "(s)");

		// Cell centers - x, y, z
		desc = "x cell centers";
		save_vector_1d(nc_file, bkg.get_x(), "x", dim2, desc, "(m)");
		desc = "y cell centers";
		save_vector_1d(nc_file, bkg.get_y(), "y", dim3, desc, "(m)");
		desc = "z cell centers";
		save_vector_1d(nc_file, bkg.get_z(), "z", dim4, desc, "(m)");
		
		// Grid edges - x, y, z
		desc = "x grid edges";
		save_vector_1d(nc_file, bkg.get_grid_x(), "grid_x", grid_dim2, desc, 
			"(m)");
		desc = "y grid edges";
		save_vector_1d(nc_file, bkg.get_grid_y(), "grid_y", grid_dim3, desc, 
			"(m)");
		desc = "z grid edges";
		save_vector_1d(nc_file, bkg.get_grid_z(), "grid_z", grid_dim4, desc, 
			"(m)");

		// Electron density
		desc = "electron density";
		save_vector_4d(nc_file, bkg.get_ne(), "electron_dens", 
			dim1, dim2, dim3, dim4, desc, "(m-3)");

		// Electron temperature
		desc = "electron temperature";
		save_vector_4d(nc_file, bkg.get_te(), "electron_temp", 
			dim1, dim2, dim3, dim4, desc, "(eV)");

		// Ion temperature
		desc = "ion temperature";
		save_vector_4d(nc_file, bkg.get_ti(), "ion_temp", 
			dim1, dim2, dim3, dim4, desc, "(eV)");

		// Plasma potential
		desc = "plasma potential";
		save_vector_4d(nc_file, bkg.get_vp(), "plasma_pot", 
			dim1, dim2, dim3, dim4, desc, "(V)");

		// Electric field
		desc = "electric field (x)";
		save_vector_4d(nc_file, bkg.get_ex(), "elec_x", 
			dim1, dim2, dim3, dim4, desc, "(V/m)");
		desc = "electric field (y)";
		save_vector_4d(nc_file, bkg.get_ey(), "elec_y", 
			dim1, dim2, dim3, dim4, desc, "(V/m)");
		desc = "electric field (z)";
		save_vector_4d(nc_file, bkg.get_ez(), "elec_z", 
			dim1, dim2, dim3, dim4, desc, "(V/m)");
			
		// Magnetic field
		desc = "magnetic field";
		save_vector_4d(nc_file, bkg.get_b(), "bmag", 
			dim1, dim2, dim3, dim4, desc, "(T)");

		// Impurity counts
		desc = "impurity counts";
		save_vector_4d(nc_file, imp_stats.get_counts(), "imp_counts", 
			dim1, dim2, dim3, dim4, desc, "");

		// Impurity density
		desc = "impurity density";
		save_vector_4d(nc_file, imp_stats.get_density(), "imp_density", 
			dim1, dim2, dim3, dim4, desc, "(m-3)");
	
		if (imp_stats.get_vel_stats())
		{
			// Impurity average x velocity
			desc = "impurity x velocity";
			save_vector_4d(nc_file, imp_stats.get_vx(), "imp_vx", 
				dim1, dim2, dim3, dim4, desc, "(m/s)");

			// Impurity average y velocity
			desc = "impurity y velocity";
			save_vector_4d(nc_file, imp_stats.get_vy(), "imp_vy", 
				dim1, dim2, dim3, dim4, desc, "(m/s)");

			// Impurity average z velocity
			desc = "impurity z velocity";
			save_vector_4d(nc_file, imp_stats.get_vz(), "imp_vz", 
				dim1, dim2, dim3, dim4, desc, "(m/s)");
		}
	}

}
