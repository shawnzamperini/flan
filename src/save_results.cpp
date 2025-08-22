/**
* @file save_results.cpp
* @brief Routines for saving results of Flan simulation
*
* Right now only outputting everything into a single netCDF file is supported.
*/
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include "background.h"
#include "flan_types.h"
#include "netcdf"
#include "options.h"
#include "save_results.h"

namespace SaveResults
{
	void save_results(const Background::Background& bkg, 
		Impurity::Statistics& imp_stats, const Options::Options& opts)
	{
		// Right now this is just a redirect to the netcdf option. I built 
		// it this way to allow other save options in the future if desired.
		save_netcdf(bkg, imp_stats, opts);
	
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
			else if constexpr (std::is_same_v<T, float>)
			{
				var = nc_file.addVar(var_name, netCDF::ncFloat, dim);
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

	// Save a 3D vector into nc_file. dims must be matched by the programmer.
	template <typename T>
	void save_vector_3d(const netCDF::NcFile& nc_file, 
		const Vectors::Vector3D<T>& vec, 
		const std::string& var_name, const netCDF::NcDim& dim1,
		const netCDF::NcDim& dim2, const netCDF::NcDim& dim3, 
		const std::string& description, const std::string& units)
	{
		if (vec.get_data().size() == 0)
		{
			std::cerr << "NetCDF Error! " << var_name << " is zero-length."
				<< " Skipping.\n";
		}
		else
		{
			// Create a 3D variable	
			//netCDF::NcVar var {nc_file.addVar(var_name, netCDF::ncDouble, 
			//		{dim1, dim2, dim3})};

			// Save as appropriate type. This isn't actually necessary since
			// netCDF will upcast a float to double, but no point in doing 
			// that if it'll just take up extra space.
			netCDF::NcVar var {};
			if constexpr (std::is_same_v<T, double>)
			{
				var = nc_file.addVar(var_name, netCDF::ncDouble, 
					{dim1, dim2, dim3});
			}
			else if constexpr (std::is_same_v<T, float>)
			{
				var = nc_file.addVar(var_name, netCDF::ncFloat, 
					{dim1, dim2, dim3});
			}
			else if constexpr (std::is_same_v<T, int>)
			{
				var = nc_file.addVar(var_name, netCDF::ncInt, 
					{dim1, dim2, dim3});
			}

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
			//netCDF::NcVar var {nc_file.addVar(var_name, netCDF::ncDouble, 
			//		{dim1, dim2, dim3, dim4})};

			// Save as appropriate type. This isn't actually necessary since
			// netCDF will upcast a float to double, but no point in doing 
			// that if it'll just take up extra space.
			netCDF::NcVar var {};
			if constexpr (std::is_same_v<T, double>)
			{
				var = nc_file.addVar(var_name, netCDF::ncDouble, 
					{dim1, dim2, dim3, dim4});
			}
			else if constexpr (std::is_same_v<T, float>)
			{
				var = nc_file.addVar(var_name, netCDF::ncFloat, 
					{dim1, dim2, dim3, dim4});
			}
			else if constexpr (std::is_same_v<T, int>)
			{
				var = nc_file.addVar(var_name, netCDF::ncInt, 
					{dim1, dim2, dim3, dim4});
			}

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

	void save_netcdf(const Background::Background& bkg, 
		Impurity::Statistics& imp_stats, const Options::Options& opts)
	{
		std::cout << "Saving netCDF file...\n";

		// Reusable string to pass in descriptions
		std::string desc {};

		// Open new netCDF file named after this case. replace tells netCDF
		// to overwrite the file if it exists.
		// NcFile's scope is private to whatever block it is defined in,
		// so we cannot pass it around functions (at least not like a normal
		// variable). The implementation details are such that it 
		// automatically closes once it goes out of scope.
		std::string nc_filename {opts.case_name() + ".nc"};

		// If file already exists, delete it otherwise we may encounter a
		// permission error if it is open elsewhere, like in a python 
		// interpretor.
		if (std::filesystem::exists(nc_filename)) 
			std::filesystem::remove(nc_filename);
		
		netCDF::NcFile nc_file {nc_filename, netCDF::NcFile::replace};
		
		// Case name
		//nc_file.putAtt("case_name", case_name);

		// Details about the simulation (system, time run, time spent, etc.)
		// To-do

		// Save all input options used for this case
		// To-do

		// Dimensions for the arrays containing data at the cell centers. 
		// First dimension is time, then x, y, z.
		netCDF::NcDim dim1 {nc_file.addDim("time", bkg.get_dim1())};
		netCDF::NcDim dim2 {nc_file.addDim("x", bkg.get_dim2())};
		netCDF::NcDim dim3 {nc_file.addDim("y", bkg.get_dim3())};
		netCDF::NcDim dim4 {nc_file.addDim("z", bkg.get_dim4())};

		// Dimensions for the grid points. They're just bigger than the
		// cell center lengths.
		netCDF::NcDim grid_dim2 {nc_file.addDim("grid x", bkg.get_dim2() + 1)};
		netCDF::NcDim grid_dim3 {nc_file.addDim("grid y", bkg.get_dim3() + 1)};
		netCDF::NcDim grid_dim4 {nc_file.addDim("grid z", bkg.get_dim4() + 1)};

		// Time
		desc = "time for each frame";
		save_vector_1d(nc_file, bkg.get_times(), "time", dim1, desc, "(s)");

		// Cell centers in computational space - x, y, z
		// We don't actually know what the dimensions will be, since it depends
		// on the computational coordinates used.
		desc = "x cell centers";
		save_vector_1d(nc_file, bkg.get_x(), "x", dim2, desc, "(?)");
		desc = "y cell centers";
		save_vector_1d(nc_file, bkg.get_y(), "y", dim3, desc, "(?)");
		desc = "z cell centers";
		save_vector_1d(nc_file, bkg.get_z(), "z", dim4, desc, "(?)");
		
		// Grid edges - x, y, z
		// We don't actually know what the dimensions will be, since it depends
		// on the computational coordinates used.
		desc = "x grid edges";
		save_vector_1d(nc_file, bkg.get_grid_x(), "grid_x", grid_dim2, desc, 
			"(?)");
		desc = "y grid edges";
		save_vector_1d(nc_file, bkg.get_grid_y(), "grid_y", grid_dim3, desc, 
			"(?)");
		desc = "z grid edges";
		save_vector_1d(nc_file, bkg.get_grid_z(), "grid_z", grid_dim4, desc, 
			"(?)");

		// Cell centers in physical space - X, Y, Z
		// Note dims are dim2, dim3 and dim4
		desc = "X cell centers";
		save_vector_3d(nc_file, bkg.get_X(), "X", dim2, dim3, dim4, desc, 
			"(m)");
		desc = "Y cell centers";
		save_vector_3d(nc_file, bkg.get_Y(), "Y", dim2, dim3, dim4, desc, 
			"(m)");
		desc = "Z cell centers";
		save_vector_3d(nc_file, bkg.get_Z(), "Z", dim2, dim3, dim4, desc, 
			"(m)");

		// Grid edges in physical space - X, Y, Z
		// Note dims are dim2, dim3 and dim4
		desc = "X grid edges";
		save_vector_3d(nc_file, bkg.get_grid_X(), "grid_X", grid_dim2, 
			grid_dim3, grid_dim4, desc, "(m)");
		desc = "Y grid edges";
		save_vector_3d(nc_file, bkg.get_grid_Y(), "grid_Y", grid_dim2, 
			grid_dim3, grid_dim4, desc, "(m)");
		desc = "Z grid edges";
		save_vector_3d(nc_file, bkg.get_grid_Z(), "grid_Z", grid_dim2, 
			grid_dim3, grid_dim4, desc, "(m)");

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
		desc = "electric field (X)";
		save_vector_4d(nc_file, bkg.get_eX(), "elec_X", 
			dim1, dim2, dim3, dim4, desc, "(V/m)");
		desc = "electric field (Y)";
		save_vector_4d(nc_file, bkg.get_eY(), "elec_Y", 
			dim1, dim2, dim3, dim4, desc, "(V/m)");
		desc = "electric field (Z)";
		save_vector_4d(nc_file, bkg.get_eZ(), "elec_Z", 
			dim1, dim2, dim3, dim4, desc, "(V/m)");

		// Background (ion) flow
		desc = "background flow (X)";
		save_vector_4d(nc_file, bkg.get_uX(), "ion_flow_X", 
			dim1, dim2, dim3, dim4, desc, "(m/s)");
		desc = "background flow (Y)";
		save_vector_4d(nc_file, bkg.get_uY(), "ion_flow_Y", 
			dim1, dim2, dim3, dim4, desc, "(m/s)");
		desc = "background flow (Z)";
		save_vector_4d(nc_file, bkg.get_uZ(), "ion_flow_Z", 
			dim1, dim2, dim3, dim4, desc, "(m/s)");
			
		// Magnetic field
		desc = "magnetic field (X)";
		save_vector_4d(nc_file, bkg.get_bX(), "bmag_X", 
			dim1, dim2, dim3, dim4, desc, "(T)");
		desc = "magnetic field (Y)";
		save_vector_4d(nc_file, bkg.get_bY(), "bmag_Y", 
			dim1, dim2, dim3, dim4, desc, "(T)");
		desc = "magnetic field (Z)";
		save_vector_4d(nc_file, bkg.get_bZ(), "bmag_Z", 
			dim1, dim2, dim3, dim4, desc, "(T)");
		desc = "magnetic field (R)";
		save_vector_4d(nc_file, bkg.get_bR(), "bmag_R", 
			dim1, dim2, dim3, dim4, desc, "(T)");

		// Magnetic field
		desc = "magnetic field gradient (X)";
		save_vector_4d(nc_file, bkg.get_gradbX(), "gradb_X", 
			dim1, dim2, dim3, dim4, desc, "(T)");
		desc = "magnetic field gradient (Y)";
		save_vector_4d(nc_file, bkg.get_gradbY(), "gradb_Y", 
			dim1, dim2, dim3, dim4, desc, "(T)");
		desc = "magnetic field gradient (Z)";
		save_vector_4d(nc_file, bkg.get_gradbZ(), "gradb_Z", 
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
			// Impurity average X velocity in physical space
			desc = "impurity X velocity";
			save_vector_4d(nc_file, imp_stats.get_vX(), "imp_vX", 
				dim1, dim2, dim3, dim4, desc, "(m/s)");

			// Impurity average Y velocity in physical space
			desc = "impurity Y velocity";
			save_vector_4d(nc_file, imp_stats.get_vY(), "imp_vY", 
				dim1, dim2, dim3, dim4, desc, "(m/s)");

			// Impurity average Z velocity in physical space
			desc = "impurity Z velocity";
			save_vector_4d(nc_file, imp_stats.get_vZ(), "imp_vZ", 
				dim1, dim2, dim3, dim4, desc, "(m/s)");
		}
		
		// Gyroradius - commenting out since it may not be useful, but
        // maybe one day it will be
		//desc = "average impurity gyroradius";
		//save_vector_4d(nc_file, imp_stats.get_gyrorad(), "imp_gyrorad", 
		//	dim1, dim2, dim3, dim4, desc, "(m)");

        // Jacobian
		desc = "Jacobian";
		save_vector_3d(nc_file, bkg.get_J(), "jacobian", dim2, 
			dim3, dim4, desc, "()");

		// Metric coefficients
		desc = "metric coefficient (00)";
		save_vector_3d(nc_file, bkg.get_gij_00(), "gij_00", dim2, 
			dim3, dim4, desc, "(?)");
		desc = "metric coefficient (01)";
		save_vector_3d(nc_file, bkg.get_gij_01(), "gij_01", dim2, 
			dim3, dim4, desc, "(?)");
		desc = "metric coefficient (02)";
		save_vector_3d(nc_file, bkg.get_gij_02(), "gij_02", dim2, 
			dim3, dim4, desc, "(?)");
		desc = "metric coefficient (11)";
		save_vector_3d(nc_file, bkg.get_gij_11(), "gij_11", dim2, 
			dim3, dim4, desc, "(?)");
		desc = "metric coefficient (12)";
		save_vector_3d(nc_file, bkg.get_gij_12(), "gij_12", dim2, 
			dim3, dim4, desc, "(?)");
		desc = "metric coefficient (22)";
		save_vector_3d(nc_file, bkg.get_gij_22(), "gij_22", dim2, 
			dim3, dim4, desc, "(?)");

		// Duals of the magnetic field
		desc = "dual basis magnetic vector (x)";
		save_vector_3d(nc_file, bkg.get_b_x(), "b_x", dim2, 
			dim3, dim4, desc, "(?)");
		desc = "dual basis magnetic vector (y)";
		save_vector_3d(nc_file, bkg.get_b_y(), "b_y", dim2, 
			dim3, dim4, desc, "(?)");
		desc = "dual basis magnetic vector (z)";
		save_vector_3d(nc_file, bkg.get_b_z(), "b_z", dim2, 
			dim3, dim4, desc, "(?)");
	}

}
