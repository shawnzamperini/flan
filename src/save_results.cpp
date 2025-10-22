/**
* @file save_results.cpp
* @brief Routines for saving results of Flan simulation
*
* Right now only outputting everything into a single netCDF file is supported.
* Important to know that NcFile inherits from NcGroup, so we can pass an
* NcFile object into a function that calls for an NcGroup.
*/
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include "background.h"
#include "config.h"
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

	// Create a NetCDF variable of the chosen type
	template <typename T>
	netCDF::NcVar create_var(const netCDF::NcGroup& nc_group,
		const std::string& var_name, const netCDF::NcDim& dim)
	{
			// Create variable with netCDF type based on which template is 
			// being used.
			netCDF::NcVar var {};
			if constexpr (std::is_same_v<T, double>)
			{
				var = nc_group.addVar(var_name, netCDF::ncDouble, dim);
			}
			else if constexpr (std::is_same_v<T, float>)
			{
				var = nc_group.addVar(var_name, netCDF::ncFloat, dim);
			}
			else if constexpr (std::is_same_v<T, int>)
			{
				var = nc_group.addVar(var_name, netCDF::ncInt, dim);
			}
			else if constexpr (std::is_same_v<T, std::vector<std::string>>)
			{
				std::cerr << "Error! If you want to save a string you must "
					<< "call save_string.\n";
			}
			else
			{
				std::cerr << "Error! Variable type for NetCDF file not "
					<< "recognized. var_name = " << var_name << '\n';
			}

			return var;
	}

	// Save a scalar to nc_group.
	template <typename T>
	void save_scalar(const netCDF::NcGroup& nc_group, const T value, 
		const std::string& var_name, const netCDF::NcDim& dim, 
		const std::string& description, const std::string& units)
	{
			// Create netCDF variable holding the templated type
			netCDF::NcVar var {create_var<T>(nc_group, var_name, dim)};

			// Put value into it (a pointer is expected here as part of the
			// netCDF API). 
			var.putVar(&value);

			// Add description
			var.putAtt("description", description);

			// Add units
			var.putAtt("units", units);
	}

	// Save a string to nc_group.
	void save_string(const netCDF::NcGroup& nc_group, const std::string value, 
		const std::string& var_name, const netCDF::NcDim& dim, 
		const std::string& description)
	{
			// Create netCDF variable holding the templated type
			netCDF::NcVar var {nc_group.addVar(var_name, netCDF::ncString, dim)};
			
			// NetCDF expects the string to be passed as an array of C-style
			// strings since the API is just a wrapper for the C API.. They just
			// can't make anything simple can they.
			std::vector<const char*> c_strs {value.c_str()};

			// Put string vector into it
			var.putVar(c_strs.data());

			// Add description
			var.putAtt("description", description);
	}

	// Save a 1D vector into nc_group. dim must be matched by the programmer.
	template <typename T>
	void save_vector_1d(const netCDF::NcGroup& nc_group, std::vector<T> vec, 
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
				var = nc_group.addVar(var_name, netCDF::ncDouble, dim);
			}
			else if constexpr (std::is_same_v<T, float>)
			{
				var = nc_group.addVar(var_name, netCDF::ncFloat, dim);
			}
			else if constexpr (std::is_same_v<T, int>)
			{
				var = nc_group.addVar(var_name, netCDF::ncInt, dim);
			}

			// Put vector data into it.
			var.putVar(vec.data());

			// Add description
			var.putAtt("description", description);

			// Add units
			var.putAtt("units", units);
		}
	}

	// Save a 3D vector into nc_group. dims must be matched by the programmer.
	template <typename T>
	void save_vector_3d(const netCDF::NcGroup& nc_group, 
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
			// Save as appropriate type. This isn't actually necessary since
			// netCDF will upcast a float to double, but no point in doing 
			// that if it'll just take up extra space.
			netCDF::NcVar var {};
			if constexpr (std::is_same_v<T, double>)
			{
				var = nc_group.addVar(var_name, netCDF::ncDouble, 
					{dim1, dim2, dim3});
			}
			else if constexpr (std::is_same_v<T, float>)
			{
				var = nc_group.addVar(var_name, netCDF::ncFloat, 
					{dim1, dim2, dim3});
			}
			else if constexpr (std::is_same_v<T, int>)
			{
				var = nc_group.addVar(var_name, netCDF::ncInt, 
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
	// Save a 4D vector into nc_group. dims must be matched by the programmer.
	template <typename T>
	void save_vector_4d(const netCDF::NcGroup& nc_group, 
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
			// Save as appropriate type. This isn't actually necessary since
			// netCDF will upcast a float to double, but no point in doing 
			// that if it'll just take up extra space.
			netCDF::NcVar var {};
			if constexpr (std::is_same_v<T, double>)
			{
				var = nc_group.addVar(var_name, netCDF::ncDouble, 
					{dim1, dim2, dim3, dim4});
			}
			else if constexpr (std::is_same_v<T, float>)
			{
				var = nc_group.addVar(var_name, netCDF::ncFloat, 
					{dim1, dim2, dim3, dim4});
			}
			else if constexpr (std::is_same_v<T, int>)
			{
				var = nc_group.addVar(var_name, netCDF::ncInt, 
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

	// Save input options to netCDF file
	void save_input_options(const netCDF::NcFile& nc_file, 
		const Options::Options& opts)
	{
		// Create group to organize things into and access dimensions
		netCDF::NcGroup input_group {nc_file.addGroup("input")};
		netCDF::NcDim dim_scalar {nc_file.getDim("scalar")};
		netCDF::NcDim dim_str {nc_file.getDim("num_string")};

		// Meta options for the simulation
		save_string(input_group, opts.case_name(), "case_name", dim_str, 
			"case name");

		// Options related to reading in Gkeyll files
		save_string(input_group, opts.gkyl_dir(), "gkyl_dir", dim_str, 
			"directory containing Gkeyll data");
		save_string(input_group, opts.gkyl_casename(), "gkyl_casename", dim_str, 
			"name of Gkeyll case");
		save_scalar(input_group, opts.gkyl_frame_start(), "gkyl_frame_start", 
			dim_scalar, "first Gkeyll frame to load", "");
		save_scalar(input_group, opts.gkyl_frame_end(), "gkyl_frame_end", 
			dim_scalar, "last Gkeyll frame to load", "");
		save_string(input_group, opts.gkyl_elec_name(), "gkyl_elec_name", 
			dim_str, "name of electron species in Gkeyll");
		save_string(input_group, opts.gkyl_ion_name(), "gkyl_ion_name", 
			dim_str, "name of deuterium species in Gkeyll");

		// Geometry options
		save_scalar(input_group, opts.lcfs_x(), "lcfs_x", 
			dim_scalar, "x value where the LCFS is located", "(?)");
		save_scalar(input_group, opts.imp_xbound_buffer(), "imp_xbound_buffer", 
			dim_scalar, "buffer for x boundary conditions", "(?)");
		save_string(input_group, opts.min_xbound_type(), "min_xbound_type", 
			dim_str, "minimum x boundary condition");

		// Impurity chracteristics
		save_scalar(input_group, opts.imp_atom_num(), "imp_atom_num", 
			dim_scalar, "impurity atomic number", "");
		save_scalar(input_group, opts.imp_mass_amu(), "imp_mass_amu", 
			dim_scalar, "impurity mass", "(amu)");
		save_scalar(input_group, opts.imp_init_charge(), "imp_init_charge", 
			dim_scalar, "impurity initial charge", "");

		// Impurity transport options
		save_scalar(input_group, opts.imp_num(), "imp_num", 
			dim_scalar, "number of primary impurities followed", "");
		save_string(input_group, opts.imp_tstart_opt(), "imp_tstart_opt", 
			dim_str, "impurity starting time option");
		save_scalar(input_group, opts.imp_tstart_val(), "imp_tstart_val", 
			dim_scalar, 
			"impurity starting time (imp_tstart_opt = single_value)", "(s)");
		save_scalar(input_group, opts.imp_trange_min(), "imp_trange_min", 
			dim_scalar, 
			"impurity starting time minimum (imp_tstart_opt = range)", "(s)");
		save_scalar(input_group, opts.imp_trange_max(), "imp_trange_max", 
			dim_scalar, 
			"impurity starting time maximum (imp_tstart_opt = range)", "(s)");
		save_string(input_group, opts.imp_xstart_opt(), "imp_xstart_opt", 
			dim_str, "impurity starting x option");
		save_scalar(input_group, opts.imp_xstart_val(), "imp_xstart_val", 
			dim_scalar, 
			"impurity starting x (imp_xstart_opt = single_value)", "(?)");
		save_scalar(input_group, opts.imp_xrange_min(), "imp_xrange_min", 
			dim_scalar, 
			"impurity starting x minimum (imp_xstart_opt = range)", "(?)");
		save_scalar(input_group, opts.imp_xrange_max(), "imp_xrange_max", 
			dim_scalar, 
			"impurity starting x maximum (imp_xstart_opt = range)", "(?)");
		save_string(input_group, opts.imp_ystart_opt(), "imp_ystart_opt", 
			dim_str, "impurity starting y option");
		save_scalar(input_group, opts.imp_ystart_val(), "imp_ystart_val", 
			dim_scalar, 
			"impurity starting y (imp_ystart_opt = single_value)", "(?)");
		save_scalar(input_group, opts.imp_yrange_min(), "imp_yrange_min", 
			dim_scalar, 
			"impurity starting y minimum (imp_ystart_opt = range)", "(?)");
		save_scalar(input_group, opts.imp_yrange_max(), "imp_yrange_max", 
			dim_scalar, 
			"impurity starting y maximum (imp_ystart_opt = range)", "(?)");
		save_string(input_group, opts.imp_zstart_opt(), "imp_zstart_opt", 
			dim_str, "impurity starting z option");
		save_scalar(input_group, opts.imp_zstart_val(), "imp_zstart_val", 
			dim_scalar, 
			"impurity starting z (imp_zstart_opt = single_value)", "(?)");
		save_scalar(input_group, opts.imp_zrange_min(), "imp_zrange_min", 
			dim_scalar, 
			"impurity starting z minimum (imp_zstart_opt = range)", "(?)");
		save_scalar(input_group, opts.imp_zrange_max(), "imp_zrange_max", 
			dim_scalar, 
			"impurity starting z maximum (imp_zstart_opt = range)", "(?)");
		save_string(input_group, opts.imp_collisions(), "imp_collisions", 
			dim_str, "impurity collisions with background plasma option");
		save_string(input_group, opts.imp_time_step_opt(), "imp_time_step_opt", 
			dim_str, "impurity time step option");
		save_scalar(input_group, opts.imp_time_step(), "imp_time_step", 
			dim_scalar, "impurity time step", "(s)");
		save_scalar(input_group, opts.imp_time_step_min(), "imp_time_step_min", 
			dim_scalar, "minimum allowed impurity time step", "(s)");
		save_scalar(input_group, opts.imp_source_scale_fact(), 
			"imp_source_scale_fact", dim_scalar, 
			"scaling factor to convert density to m-3", "(particles/s)");
		save_string(input_group, opts.imp_iz_recomb(), "imp_iz_recomb", 
			dim_str, "impurity ionization/recombination option");

		// Variance reduction options
		save_string(input_group, opts.var_red_split(), "var_red_split", 
			dim_str, "variance reduction particle splitting option");
		save_string(input_group, opts.var_red_import(), "var_red_import", 
			dim_str, "importance criteria for variance reduction");
		save_scalar(input_group, opts.var_red_freq(), "var_red_freq", 
			dim_scalar, "variance reduction criteria update frequency", "");
		save_scalar(input_group, opts.var_red_min_weight(), 
			"var_red_min_weight", dim_scalar, 
			"minimum allowed particle weight for variance reduction", "");
		save_scalar(input_group, opts.var_red_med_mod(), 
			"var_red_med_mod", dim_scalar, 
			"multiplier for median importance (var_red_import = median)", "");
		save_string(input_group, opts.var_red_rusrol(), "var_red_rusrol", 
			dim_str, "Russian roullette variance reduciton option");
		save_scalar(input_group, opts.var_red_rusrol_prob(), 
			"var_red_rusrol_prob", dim_scalar, 
			"kill probability for Russian roullette variance reduction", "");

		// OpenADAS options
		save_string(input_group, opts.openadas_root(), "openadas_root", 
			dim_scalar, "directory containing OpenADAS data");
		save_scalar(input_group, opts.openadas_year(), "openadas_year", 
			dim_scalar, "year for OpenADAS files", "");

	}

	void save_geometry(const netCDF::NcFile& nc_file, 
		const Background::Background& bkg)
	{
		// Create group to organize things into and access dimensions
		netCDF::NcGroup geo_group {nc_file.addGroup("geometry")};
		netCDF::NcDim dim1 {nc_file.getDim("time")};
		netCDF::NcDim dim2 {nc_file.getDim("x")};
		netCDF::NcDim dim3 {nc_file.getDim("y")};
		netCDF::NcDim dim4 {nc_file.getDim("z")};
		netCDF::NcDim grid_dim2 {nc_file.getDim("grid x")};
		netCDF::NcDim grid_dim3 {nc_file.getDim("grid y")};
		netCDF::NcDim grid_dim4 {nc_file.getDim("grid z")};

		// Time
		save_vector_1d(geo_group, bkg.get_times(), "time", dim1, 
			"time for each frame", "(s)");

		// Cell centers in computational space - x, y, z
		// We don't actually know what the units will be, since it depends
		// on the computational coordinates used.
		save_vector_1d(geo_group, bkg.get_x(), "x", dim2, "x cell centers", 
			"(?)");
		save_vector_1d(geo_group, bkg.get_y(), "y", dim3, "y cell centers", 
			"(?)");
		save_vector_1d(geo_group, bkg.get_z(), "z", dim4, "z cell centers", 
			"(?)");

		// Grid edges - x, y, z
		// We don't actually know what the dimensions will be, since it depends
		// on the computational coordinates used.
		save_vector_1d(geo_group, bkg.get_grid_x(), "grid_x", grid_dim2, 
			"x grid edges", "(?)");
		save_vector_1d(geo_group, bkg.get_grid_y(), "grid_y", grid_dim3, 
			"y grid edges", "(?)");
		save_vector_1d(geo_group, bkg.get_grid_z(), "grid_z", grid_dim4, 
			"z grid edges", "(?)");

		// Cell centers in physical space - X, Y, Z
		// Note dims are dim2, dim3 and dim4
		save_vector_3d(geo_group, bkg.get_X(), "X", dim2, dim3, dim4, 
			"X cell centers", "(m)");
		save_vector_3d(geo_group, bkg.get_Y(), "Y", dim2, dim3, dim4, 
			"Y cell centers", "(m)");
		save_vector_3d(geo_group, bkg.get_Z(), "Z", dim2, dim3, dim4, 
			"Z cell centers", "(m)");

		// Grid edges in physical space - X, Y, Z
		// Note dims are dim2, dim3 and dim4
		save_vector_3d(geo_group, bkg.get_grid_X(), "grid_X", grid_dim2, 
			grid_dim3, grid_dim4, "X grid edges", "(m)");
		save_vector_3d(geo_group, bkg.get_grid_Y(), "grid_Y", grid_dim2, 
			grid_dim3, grid_dim4, "Y grid edges", "(m)");
		save_vector_3d(geo_group, bkg.get_grid_Z(), "grid_Z", grid_dim2, 
			grid_dim3, grid_dim4, "Z grid edges", "(m)");

		// Jacobian
		save_vector_3d(geo_group, bkg.get_J(), "J", dim2, 
			dim3, dim4, "Jacobian", "()");

		// Metric coefficients
		save_vector_3d(geo_group, bkg.get_gij_00(), "gij_00", dim2, 
			dim3, dim4, "metric coefficient (00)", "(?)");
		save_vector_3d(geo_group, bkg.get_gij_01(), "gij_01", dim2, 
			dim3, dim4, "metric coefficient (01)", "(?)");
		save_vector_3d(geo_group, bkg.get_gij_02(), "gij_02", dim2, 
			dim3, dim4, "metric coefficient (02)", "(?)");
		save_vector_3d(geo_group, bkg.get_gij_11(), "gij_11", dim2, 
			dim3, dim4, "metric coefficient (11)", "(?)");
		save_vector_3d(geo_group, bkg.get_gij_12(), "gij_12", dim2, 
			dim3, dim4, "metric coefficient (12)", "(?)");
		save_vector_3d(geo_group, bkg.get_gij_22(), "gij_22", dim2, 
			dim3, dim4, "metric coefficient (22)", "(?)");

		// Duals of the magnetic field
		save_vector_3d(geo_group, bkg.get_b_x(), "b_x", dim2, 
			dim3, dim4, "dual basis magnetic vector (x)", "(?)");
		save_vector_3d(geo_group, bkg.get_b_y(), "b_y", dim2, 
			dim3, dim4, "dual basis magnetic vector (y)", "(?)");
		save_vector_3d(geo_group, bkg.get_b_z(), "b_z", dim2, 
			dim3, dim4, "dual basis magnetic vector (z)", "(?)");
	}

	void save_background(const netCDF::NcFile& nc_file, 
		const Options::Options& opts, const Background::Background& bkg)
	{
		// Create group to organize things into and access dimensions
		netCDF::NcGroup bkg_group {nc_file.addGroup("background")};
		netCDF::NcDim dim1 {nc_file.getDim("time")};
		netCDF::NcDim dim2 {nc_file.getDim("x")};
		netCDF::NcDim dim3 {nc_file.getDim("y")};
		netCDF::NcDim dim4 {nc_file.getDim("z")};

		// Electron density
		save_vector_4d(bkg_group, bkg.get_ne(), "ne", 
			dim1, dim2, dim3, dim4, "electron density", "(m-3)");

		// Electron temperature
		save_vector_4d(bkg_group, bkg.get_te(), "Te", 
			dim1, dim2, dim3, dim4, "electron temperature", "(eV)");

		// Ion temperature
		save_vector_4d(bkg_group, bkg.get_ti(), "Ti", 
			dim1, dim2, dim3, dim4, "ion temperature", "(eV)");

		// Plasma potential
		save_vector_4d(bkg_group, bkg.get_vp(), "Vp", 
			dim1, dim2, dim3, dim4, "plasma potential", "(V)");

		// Electric field
		save_vector_4d(bkg_group, bkg.get_eX(), "E_X", 
			dim1, dim2, dim3, dim4, "electric field (X)", "(V/m)");
		save_vector_4d(bkg_group, bkg.get_eY(), "E_Y", 
			dim1, dim2, dim3, dim4, "electric field (Y)", "(V/m)");
		save_vector_4d(bkg_group, bkg.get_eZ(), "E_Z", 
			dim1, dim2, dim3, dim4, "electric field (Z)", "(V/m)");

		// Background (ion) flow
		save_vector_4d(bkg_group, bkg.get_uX(), "Ui_X", 
			dim1, dim2, dim3, dim4, "background flow (X)", "(m/s)");
		save_vector_4d(bkg_group, bkg.get_uY(), "Ui_Y", 
			dim1, dim2, dim3, dim4, "background flow (Y)", "(m/s)");
		save_vector_4d(bkg_group, bkg.get_uZ(), "Ui_Z", 
			dim1, dim2, dim3, dim4, "background flow (Z)", "(m/s)");

		// Magnetic field
		save_vector_4d(bkg_group, bkg.get_bX(), "B_X", 
			dim1, dim2, dim3, dim4, "magnetic field (X)", "(T)");
		save_vector_4d(bkg_group, bkg.get_bY(), "B_Y", 
			dim1, dim2, dim3, dim4, "magnetic field (Y)", "(T)");
		save_vector_4d(bkg_group, bkg.get_bZ(), "B_Z", 
			dim1, dim2, dim3, dim4, "magnetic field (Z)", "(T)");
		save_vector_4d(bkg_group, bkg.get_bR(), "B_R", 
			dim1, dim2, dim3, dim4, "magnetic field (R)", "(T)");

		// Magnetic field gradient
		save_vector_4d(bkg_group, bkg.get_gradbX(), "gradB_X", 
			dim1, dim2, dim3, dim4, "magnetic field gradient (X)", "(T)");
		save_vector_4d(bkg_group, bkg.get_gradbY(), "gradB_Y", 
			dim1, dim2, dim3, dim4, "magnetic field gradient (Y)", "(T)");
		save_vector_4d(bkg_group, bkg.get_gradbZ(), "gradB_Z", 
			dim1, dim2, dim3, dim4, "magnetic field gradient (Z)", "(T)");

		// Electric field gradient
		if (opts.calc_grad_elec_int() == 1)
		{
			save_vector_4d(bkg_group, bkg.get_gradeX(), "gradE_X", dim1, dim2, 
				dim3, dim4, "electric field gradient (X)", "(V/m2)");
			save_vector_4d(bkg_group, bkg.get_gradeY(), "gradE_Y", dim1, dim2, 
				dim3, dim4, "electric field gradient (Y)", "(V/m2)");
			save_vector_4d(bkg_group, bkg.get_gradeZ(), "gradE_Z", dim1, dim2, 
				dim3, dim4, "electric field gradient (Z)", "(V/m2)");
		}

	}

	void save_output(const netCDF::NcFile& nc_file, 
		const Options::Options& opts, const Background::Background& bkg,
		Impurity::Statistics& imp_stats)
	{
		// Create group to organize things into and access dimensions
		netCDF::NcGroup out_group {nc_file.addGroup("output")};
		netCDF::NcDim dim1 {nc_file.getDim("time")};
		netCDF::NcDim dim2 {nc_file.getDim("x")};
		netCDF::NcDim dim3 {nc_file.getDim("y")};
		netCDF::NcDim dim4 {nc_file.getDim("z")};

		// Impurity counts
		save_vector_4d(out_group, imp_stats.get_counts(), "Nz", 
			dim1, dim2, dim3, dim4, "impurity counts", "");

		// Impurity density
		save_vector_4d(out_group, imp_stats.get_density(), "nz", 
			dim1, dim2, dim3, dim4, "impurity density", "(m-3)");

		// Impurity average X velocity in physical space
		save_vector_4d(out_group, imp_stats.get_vX(), "v_X", 
			dim1, dim2, dim3, dim4, "average impurity X velocity", "(m/s)");

		// Impurity average Y velocity in physical space
		save_vector_4d(out_group, imp_stats.get_vY(), "v_Y", 
			dim1, dim2, dim3, dim4, "average impurity Y velocity", "(m/s)");

		// Impurity average Z velocity in physical space
		save_vector_4d(out_group, imp_stats.get_vZ(), "v_Z", 
			dim1, dim2, dim3, dim4, "average impurity Z velocity", "(m/s)");

		// Impurity charge
		save_vector_4d(out_group, imp_stats.get_charge(), "qz", 
			dim1, dim2, dim3, dim4, "average impurity charge", "()");
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

		// Scalars just have length 1 dimension.
		netCDF::NcDim dim_scalar {nc_file.addDim("scalar", 1)};

		// Dimension for string variables. Also length one (technically can
		// reuse dim_scalar, but leave separate in case this changes).
		netCDF::NcDim dim_str {nc_file.addDim("num_string", 1)};

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

		// Details about the simulation (system, time run, time spent, etc.)
		// To-do

		// Save commit hash.
		const std::string commit_hash {GIT_COMMIT_HASH};
		desc = "commit hash";
		save_string(nc_file, commit_hash, "commit", dim_str, desc);

		// Save all input options used for this case
		save_input_options(nc_file, opts);

		// Save all geometry related variables
		save_geometry(nc_file, bkg);

		// Save all background related variables
		save_background(nc_file, opts, bkg);

		// Save impurity transport results
		save_output(nc_file, opts, bkg, imp_stats);

		std::cout << "Saved as " << nc_filename << '\n';
	}

}
