/**
* @file save_results.cpp
* @brief Routines for saving results of Flan simulation
*
* Right now only outputting everything into a single netCDF file is supported.
* Important to know that NcFile inherits from NcGroup, so we can pass an
* NcFile object into a function that calls for an NcGroup.
*
* The file is named <case_name>.nc and has four groups: input (the options
* used), geometry (grids, metric coefficients and basis vectors),
* background (the plasma) and output (the impurity results). The 4D arrays
* have dimensions (time, x, y, z), and the grid edge arrays use the
* "grid x", "grid y" and "grid z" dimensions, one longer than x, y and z.
* Every variable has a description attribute, and numeric ones a units
* attribute. Units of "(?)" depend on the computational coordinates used.
*
* A background with more than one zone (e.g., SOLEDGE3X) is saved
* differently, since each zone has its own grid. The input group and the
* times (geometry/time) are the same. A zones group holds the number of
* cells in each zone and what is across each zone face. Then each zone has
* its own group, zone1, zone2, ... (numbered from 1), with its own x, y, z
* and grid dimensions and geometry, background and output subgroups. These
* hold the same variables as for a single zone, interior cells only (no
* ghost cells), plus the wall mask chi and the poloidal cell corners
* R_corners and Z_corners for plotting.
*/
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include "background.h"
#include "config.h"
#include "errors.h"
#include "flan_types.h"
#include "netcdf"
#include "options.h"
#include "save_results.h"
#include "vectors.h"
#include "zones.h"

namespace SaveResults
{
	/**
	* Saves the results of a simulation.
	*
	* Entry point for saving. Right now this just calls save_netcdf, but it
	* leaves room for other output formats.
	*
	* @param bkg Background plasma, saved along with the grids and geometry.
	* @param imp_stats Impurity statistics to save.
	* @param opts Simulation options, saved in the input group.
	*/
	void save_results(const Background::Background& bkg, 
		Impurity::Statistics& imp_stats, const Options::Options& opts)
	{
		// Right now this is just a redirect to the netcdf option. I built 
		// it this way to allow other save options in the future if desired.
		save_netcdf(bkg, imp_stats, opts);
	
	}

	/**
	* Creates a netCDF variable with the netCDF type matching T.
	*
	* double, float and int are supported. Strings must be saved with
	* save_string instead. Any other type is a programming error and stops
	* the run with Errors::fatal_error.
	*
	* @tparam T C++ type of the values that will be stored.
	* @param nc_group Group (or file) to add the variable to.
	* @param var_name Name of the variable.
	* @param dim Dimension of the variable.
	* @return The new, empty variable.
	*/
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
				Errors::fatal_error("If you want to save a string you must "
					"call save_string. var_name = " + var_name);
			}
			else
			{
				Errors::fatal_error("Variable type for NetCDF file not "
					"recognized. var_name = " + var_name);
			}

			return var;
	}

	/**
	* Saves a single number as a variable in nc_group.
	*
	* @tparam T Type of the value: double, float or int.
	* @param nc_group Group (or file) to save into.
	* @param value Value to save.
	* @param var_name Name of the variable.
	* @param dim Length 1 dimension for scalars ("scalar").
	* @param description Description attribute.
	* @param units Units attribute, e.g., "(s)". Empty for unitless values.
	*/
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

	/**
	* Saves a string as a variable in nc_group.
	*
	* @param nc_group Group (or file) to save into.
	* @param value String to save.
	* @param var_name Name of the variable.
	* @param dim Length 1 dimension for strings ("num_string").
	* @param description Description attribute.
	*/
	void save_string(const netCDF::NcGroup& nc_group, const std::string value, 
		const std::string& var_name, const netCDF::NcDim& dim, 
		const std::string& description)
	{
			// Create netCDF string variable
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

	/**
	* Saves a 1D vector as a variable in nc_group.
	*
	* An empty vector is skipped with a warning, since netCDF can't store a
	* zero-length variable. dim must have the same length as vec.
	*
	* @tparam T Type of the values: double, float or int.
	* @param nc_group Group (or file) to save into.
	* @param vec Values to save.
	* @param var_name Name of the variable.
	* @param dim Dimension of the variable, the same length as vec.
	* @param description Description attribute.
	* @param units Units attribute.
	*/
	template <typename T>
	void save_vector_1d(const netCDF::NcGroup& nc_group, std::vector<T> vec, 
		const std::string& var_name, const netCDF::NcDim& dim, 
		const std::string& description, const std::string& units)
	{
		// Vector cannot be zero-length.
		if (vec.size() == 0)
		{
			std::cerr << "Warning! " << var_name << " is zero-length."
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

	/**
	* Saves a Vector3D, or a single-zone ZonedVector without ghost cells, as
	* a variable in nc_group.
	*
	* The flattened data is written as is, since netCDF uses the same
	* row-major order. An empty vector is skipped with a warning. The
	* dimensions must match the shape of vec.
	*
	* @tparam V Type of the vector. Its values must be double, float or
	*         int.
	* @param nc_group Group (or file) to save into.
	* @param vec Values to save.
	* @param var_name Name of the variable.
	* @param dim1 First dimension of the variable.
	* @param dim2 Second dimension of the variable.
	* @param dim3 Third dimension of the variable.
	* @param description Description attribute.
	* @param units Units attribute.
	*/
	template <typename V>
	void save_vector_3d(const netCDF::NcGroup& nc_group, const V& vec,
		const std::string& var_name, const netCDF::NcDim& dim1,
		const netCDF::NcDim& dim2, const netCDF::NcDim& dim3, 
		const std::string& description, const std::string& units)
	{
		// Type of the values
		using T = std::remove_cvref_t<decltype(vec.get_data()[0])>;

		if (vec.get_data().size() == 0)
		{
			std::cerr << "Warning! " << var_name << " is zero-length."
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

	/**
	* Saves a Vector4D, or a single-zone ZonedVector without ghost cells, as
	* a variable in nc_group.
	*
	* The flattened data is written as is, since netCDF uses the same
	* row-major order. An empty vector is skipped with a warning. The
	* dimensions must match the shape of vec.
	*
	* @tparam V Type of the vector. Its values must be double, float or
	*         int.
	* @param nc_group Group (or file) to save into.
	* @param vec Values to save.
	* @param var_name Name of the variable.
	* @param dim1 First dimension of the variable (usually time).
	* @param dim2 Second dimension of the variable.
	* @param dim3 Third dimension of the variable.
	* @param dim4 Fourth dimension of the variable.
	* @param description Description attribute.
	* @param units Units attribute.
	*/
	template <typename V>
	void save_vector_4d(const netCDF::NcGroup& nc_group, const V& vec,
		const std::string& var_name, const netCDF::NcDim& dim1,
		const netCDF::NcDim& dim2, const netCDF::NcDim& dim3, 
		const netCDF::NcDim& dim4, const std::string& description, 
		const std::string& units)
	{
		// Type of the values
		using T = std::remove_cvref_t<decltype(vec.get_data()[0])>;

		if (vec.get_data().size() == 0)
		{
			std::cerr << "Warning! " << var_name << " is zero-length."
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

	/**
	* Saves every input option into the "input" group.
	*
	* Each option is saved under its option name, with a description and,
	* for numbers, units. The groups follow the Input Options page of the
	* docs. An option added to load_input should be added here too.
	*
	* @param nc_file File to save into. The "scalar" and "num_string"
	*        dimensions must already exist.
	* @param opts Simulation options.
	*/
	void save_input_options(const netCDF::NcFile& nc_file, 
		const Options::Options& opts)
	{
		// Create group to organize things into and access dimensions
		netCDF::NcGroup input_group {nc_file.addGroup("input")};
		netCDF::NcDim dim_scalar {nc_file.getDim("scalar")};
		netCDF::NcDim dim_str {nc_file.getDim("num_string")};

		// Simulation
		save_string(input_group, opts.case_name(), "case_name", dim_str, 
			"case name");
		save_string(input_group, opts.use_gpu(), "use_gpu", dim_str, 
			"GPU acceleration option");
		save_scalar(input_group, opts.seed(), "seed", dim_scalar, 
			"base seed for the random number generators (negative = from "
			"the clock)", "");
		save_scalar(input_group, opts.slot_cap(), "slot_cap", dim_scalar, 
			"maximum number of particle slots", "");
		save_scalar(input_group, opts.print_interval(), "print_interval", 
			dim_scalar, "number of progress updates printed", "");

		// Background plasma
		save_string(input_group, opts.bkg_source(), "bkg_source", dim_str, 
			"source of background plasma");
		save_scalar(input_group, opts.min_ne(), "min_ne", dim_scalar, 
			"lower limit on the electron density", "(m-3)");
		save_scalar(input_group, opts.min_te(), "min_te", dim_scalar, 
			"lower limit on the electron temperature", "(eV)");
		save_scalar(input_group, opts.min_ti(), "min_ti", dim_scalar, 
			"lower limit on the ion temperature", "(eV)");

		// Testing
		save_string(input_group, opts.test_opt(), "test_opt", dim_str, 
			"test background plasma option");
		save_string(input_group, opts.save_track(), "save_track", dim_str, 
			"save particle t,x,y,z track option");

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
			dim_str, "name of main ion species in Gkeyll");
		save_scalar(input_group, opts.gkyl_elec_mass_amu(), 
			"gkyl_elec_mass_amu", dim_scalar, 
			"electron mass used in Gkeyll", "(amu)");
		save_scalar(input_group, opts.gkyl_ion_mass_amu(), 
			"gkyl_ion_mass_amu", dim_scalar, 
			"main ion mass used in Gkeyll", "(amu)");
		save_string(input_group, opts.gkyl_file_type(), "gkyl_file_type", 
			dim_str, "file format of the Gkeyll output");
		save_string(input_group, opts.gkyl_moment_type(), "gkyl_moment_type", 
			dim_str, "type of moment files in the Gkeyll output");
		save_string(input_group, opts.calc_grad_elec(), "calc_grad_elec", 
			dim_str, "calculate the electric field gradient option");

		// Options related to reading in SOLEDGE3X files
		save_string(input_group, opts.s3x_dir(), "s3x_dir", dim_str,
			"directory containing SOLEDGE3X data");
		save_string(input_group, opts.s3x_mesh_file(), "s3x_mesh_file", dim_str,
			"name of SOLEDGE3X mesh file");
		save_scalar(input_group, opts.s3x_frame_start(), "s3x_frame_start", dim_scalar,
			"first SOLEDGE3X plasma frame to load", "");
		save_scalar(input_group, opts.s3x_frame_end(), "s3x_frame_end", dim_scalar,
			"last SOLEDGE3X plasma frame to load", "");
		save_string(input_group, opts.s3x_ref_file(), "s3x_ref_file", dim_str,
			"name of SOLEDGE3X reference (normalization) file");
		save_string(input_group, opts.s3x_chi_absorb(), "s3x_chi_absorb", dim_str,
			"absorb particles inside the SOLEDGE3X wall option");

		// Geometry and boundaries
		save_string(input_group, opts.tbound_type(), "tbound_type", 
			dim_str, "time boundary condition");
		save_string(input_group, opts.min_xbound_type(), "min_xbound_type", 
			dim_str, "minimum x boundary condition");
		save_string(input_group, opts.max_xbound_type(), "max_xbound_type", 
			dim_str, "maximum x boundary condition");
		save_string(input_group, opts.min_ybound_type(), "min_ybound_type", 
			dim_str, "minimum y boundary condition");
		save_string(input_group, opts.max_ybound_type(), "max_ybound_type", 
			dim_str, "maximum y boundary condition");
		save_string(input_group, opts.min_zbound_type(), "min_zbound_type", 
			dim_str, "minimum z boundary condition");
		save_string(input_group, opts.max_zbound_type(), "max_zbound_type", 
			dim_str, "maximum z boundary condition");
		save_scalar(input_group, opts.imp_xbound_buffer(), "imp_xbound_buffer", 
			dim_scalar, "buffer for x boundary conditions", "(?)");
		save_scalar(input_group, opts.imp_ybound_buffer(), "imp_ybound_buffer", 
			dim_scalar, "buffer for y boundary conditions", "(?)");
		save_scalar(input_group, opts.imp_zbound_buffer(), "imp_zbound_buffer", 
			dim_scalar, "buffer for z boundary conditions", "(?)");
		save_scalar(input_group, opts.lcfs_x(), "lcfs_x", 
			dim_scalar, "x value where the LCFS is located", "(?)");
		save_scalar(input_group, opts.lcfs_y(), "lcfs_y", 
			dim_scalar, "y value where the LCFS is located", "(?)");
		save_scalar(input_group, opts.lcfs_z(), "lcfs_z", 
			dim_scalar, "z value where the LCFS is located", "(?)");
		save_scalar(input_group, opts.sep_x_bc_xp_z1(), "sep_x_bc_xp_z1", 
			dim_scalar, 
			"first z value for the separatrix x boundary near the X-point",
			"(?)");
		save_scalar(input_group, opts.sep_x_bc_xp_z2(), "sep_x_bc_xp_z2", 
			dim_scalar, 
			"second z value for the separatrix x boundary near the X-point",
			"(?)");

		// Impurity characteristics
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
		save_scalar(input_group, opts.imp_zone_start_val(),
			"imp_zone_start_val", dim_scalar,
			"zone impurities start in (numbered from 1)", "");
		save_scalar(input_group, opts.imp_zrange_min(), "imp_zrange_min", 
			dim_scalar, 
			"impurity starting z minimum (imp_zstart_opt = range)", "(?)");
		save_scalar(input_group, opts.imp_zrange_max(), "imp_zrange_max", 
			dim_scalar, 
			"impurity starting z maximum (imp_zstart_opt = range)", "(?)");
		save_string(input_group, opts.imp_temp_start_opt(), 
			"imp_temp_start_opt", dim_str, 
			"impurity starting temperature option");
		save_scalar(input_group, opts.imp_temp_start_val(), 
			"imp_temp_start_val", dim_scalar, 
			"impurity starting temperature (imp_temp_start_opt = "
			"single_value)", "(eV)");
		save_string(input_group, opts.imp_collisions(), "imp_collisions", 
			dim_str, "impurity collisions with background plasma option");
		save_string(input_group, opts.imp_time_step_opt(), "imp_time_step_opt", 
			dim_str, "impurity time step option");
		save_scalar(input_group, opts.imp_time_step(), "imp_time_step", 
			dim_scalar, "impurity time step", "(s)");
		save_scalar(input_group, opts.imp_time_step_min(), "imp_time_step_min", 
			dim_scalar, "minimum allowed impurity time step", "(s)");
		save_scalar(input_group, opts.imp_max_time(), "imp_max_time", 
			dim_scalar, "maximum time a particle is followed", "(s)");
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
			dim_str, "Russian roulette variance reduction option");
		save_scalar(input_group, opts.var_red_rusrol_prob(), 
			"var_red_rusrol_prob", dim_scalar, 
			"kill probability for Russian roulette variance reduction", "");

		// OpenADAS options
		save_string(input_group, opts.openadas_root(), "openadas_root", 
			dim_str, "directory containing OpenADAS data");
		save_scalar(input_group, opts.openadas_year(), "openadas_year", 
			dim_scalar, "year for OpenADAS files", "");

	}

	/**
	* Saves the grids and geometry into the "geometry" group.
	*
	* This is the background times, the computational cell centers and
	* edges, the physical (X, Y, Z) cell centers and edges, the Jacobian,
	* the metric coefficients, the tangent and reciprocal basis vectors, and
	* the dual basis components of the magnetic field.
	*
	* @param nc_file File to save into. The time, x, y, z and grid
	*        dimensions must already exist.
	* @param bkg Background plasma holding the geometry.
	*/
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
		save_vector_1d(geo_group, bkg.get_x(0), "x", dim2, "x cell centers", 
			"(?)");
		save_vector_1d(geo_group, bkg.get_y(0), "y", dim3, "y cell centers", 
			"(?)");
		save_vector_1d(geo_group, bkg.get_z(0), "z", dim4, "z cell centers", 
			"(?)");

		// Grid edges - x, y, z
		// We don't actually know what the units will be, since it depends
		// on the computational coordinates used.
		save_vector_1d(geo_group, bkg.get_grid_x(0), "grid_x", grid_dim2, 
			"x grid edges", "(?)");
		save_vector_1d(geo_group, bkg.get_grid_y(0), "grid_y", grid_dim3, 
			"y grid edges", "(?)");
		save_vector_1d(geo_group, bkg.get_grid_z(0), "grid_z", grid_dim4, 
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
		// Note dims are the grid dims, one longer than dim2, dim3 and dim4
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

		// Tangent basis vector components
		save_vector_3d(geo_group, bkg.get_dXdx(), "dXdx", dim2, 
			dim3, dim4, "tangent basis vector (dX/dx)", "(?)");
		save_vector_3d(geo_group, bkg.get_dYdx(), "dYdx", dim2, 
			dim3, dim4, "tangent basis vector (dY/dx)", "(?)");
		save_vector_3d(geo_group, bkg.get_dZdx(), "dZdx", dim2, 
			dim3, dim4, "tangent basis vector (dZ/dx)", "(?)");
		save_vector_3d(geo_group, bkg.get_dXdy(), "dXdy", dim2, 
			dim3, dim4, "tangent basis vector (dX/dy)", "(?)");
		save_vector_3d(geo_group, bkg.get_dYdy(), "dYdy", dim2, 
			dim3, dim4, "tangent basis vector (dY/dy)", "(?)");
		save_vector_3d(geo_group, bkg.get_dZdy(), "dZdy", dim2, 
			dim3, dim4, "tangent basis vector (dZ/dy)", "(?)");
		save_vector_3d(geo_group, bkg.get_dXdz(), "dXdz", dim2, 
			dim3, dim4, "tangent basis vector (dX/dz)", "(?)");
		save_vector_3d(geo_group, bkg.get_dYdz(), "dYdz", dim2, 
			dim3, dim4, "tangent basis vector (dY/dz)", "(?)");
		save_vector_3d(geo_group, bkg.get_dZdz(), "dZdz", dim2, 
			dim3, dim4, "tangent basis vector (dZ/dz)", "(?)");

		// Reciprocal basis vector components
		save_vector_3d(geo_group, bkg.get_dxdX(), "dxdX", dim2, 
			dim3, dim4, "reciprocal basis vector (dx/dX)", "(?)");
		save_vector_3d(geo_group, bkg.get_dxdY(), "dxdY", dim2, 
			dim3, dim4, "reciprocal basis vector (dx/dY)", "(?)");
		save_vector_3d(geo_group, bkg.get_dxdZ(), "dxdZ", dim2, 
			dim3, dim4, "reciprocal basis vector (dx/dZ)", "(?)");
		save_vector_3d(geo_group, bkg.get_dydX(), "dydX", dim2, 
			dim3, dim4, "reciprocal basis vector (dy/dX)", "(?)");
		save_vector_3d(geo_group, bkg.get_dydY(), "dydY", dim2, 
			dim3, dim4, "reciprocal basis vector (dy/dY)", "(?)");
		save_vector_3d(geo_group, bkg.get_dydZ(), "dydZ", dim2, 
			dim3, dim4, "reciprocal basis vector (dy/dZ)", "(?)");
		save_vector_3d(geo_group, bkg.get_dzdX(), "dzdX", dim2, 
			dim3, dim4, "reciprocal basis vector (dz/dX)", "(?)");
		save_vector_3d(geo_group, bkg.get_dzdY(), "dzdY", dim2, 
			dim3, dim4, "reciprocal basis vector (dz/dY)", "(?)");
		save_vector_3d(geo_group, bkg.get_dzdZ(), "dzdZ", dim2, 
			dim3, dim4, "reciprocal basis vector (dz/dZ)", "(?)");

		// Duals of the magnetic field
		save_vector_3d(geo_group, bkg.get_b_x(), "b_x", dim2, 
			dim3, dim4, "dual basis magnetic vector (x)", "(?)");
		save_vector_3d(geo_group, bkg.get_b_y(), "b_y", dim2, 
			dim3, dim4, "dual basis magnetic vector (y)", "(?)");
		save_vector_3d(geo_group, bkg.get_b_z(), "b_z", dim2, 
			dim3, dim4, "dual basis magnetic vector (z)", "(?)");
	}

	/**
	* Saves the background plasma into the "background" group.
	*
	* This is the density, temperatures, plasma potential, electric field,
	* ion flow, magnetic field and its gradient, all with dimensions (time,
	* x, y, z). The electric field gradient is only saved if it was
	* calculated (calc_grad_elec).
	*
	* @param nc_file File to save into. The time, x, y and z dimensions
	*        must already exist.
	* @param opts Simulation options, used for calc_grad_elec.
	* @param bkg Background plasma to save.
	*/
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
		save_vector_4d(bkg_group, bkg.get_gradbX(), "dBdX", 
			dim1, dim2, dim3, dim4, "magnetic field gradient (X)", "(T/m)");
		save_vector_4d(bkg_group, bkg.get_gradbY(), "dBdY", 
			dim1, dim2, dim3, dim4, "magnetic field gradient (Y)", "(T/m)");
		save_vector_4d(bkg_group, bkg.get_gradbZ(), "dBdZ", 
			dim1, dim2, dim3, dim4, "magnetic field gradient (Z)", "(T/m)");

		// Electric field gradient
		if (opts.calc_grad_elec_int() == 1)
		{
			save_vector_4d(bkg_group, bkg.get_gradeX(), "dEdX", dim1, dim2, 
				dim3, dim4, "electric field gradient (X)", "(V/m2)");
			save_vector_4d(bkg_group, bkg.get_gradeY(), "dEdY", dim1, dim2, 
				dim3, dim4, "electric field gradient (Y)", "(V/m2)");
			save_vector_4d(bkg_group, bkg.get_gradeZ(), "dEdZ", dim1, dim2, 
				dim3, dim4, "electric field gradient (Z)", "(V/m2)");
		}

	}

	/**
	* Saves the impurity results into the "output" group.
	*
	* This is the impurity counts and density, the average Cartesian and
	* computational velocities, the average charge and the average Nanbu
	* collisionality s, all with dimensions (time, x, y, z). The particle
	* track is saved too if save_track is on.
	*
	* @param nc_file File to save into. The time, x, y, z and track
	*        dimensions must already exist.
	* @param opts Simulation options, used for save_track.
	* @param bkg Background plasma. Currently unused.
	* @param imp_stats Impurity statistics to save.
	*/
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

		// Impurity average x velocity in computational space
		save_vector_4d(out_group, imp_stats.get_vx(), "v_x", 
			dim1, dim2, dim3, dim4, "average impurity x velocity", "(?/s)");
			
		// Impurity average y velocity in computational space
		save_vector_4d(out_group, imp_stats.get_vy(), "v_y", 
			dim1, dim2, dim3, dim4, "average impurity y velocity", "(?/s)");

		// Impurity average z velocity in computational space
		save_vector_4d(out_group, imp_stats.get_vz(), "v_z", 
			dim1, dim2, dim3, dim4, "average impurity z velocity", "(?/s)");

		// Impurity charge
		save_vector_4d(out_group, imp_stats.get_charge(), "qz", 
			dim1, dim2, dim3, dim4, "average impurity charge", "()");

		// Nanbu impurity-ion collisionality strength factor
		save_vector_4d(out_group, imp_stats.get_s(), "nanbu_s", 
			dim1, dim2, dim3, dim4, "average Nanbu imp-ion collisionality strength", 
			"()");

		// Impurity track
		if (opts.save_track_int() == 1)
		{
			netCDF::NcDim track_dim {nc_file.getDim("track")};
			save_vector_1d(out_group, imp_stats.get_track_t(), "track_t", 
				track_dim, "t particle track", "(?)");
			save_vector_1d(out_group, imp_stats.get_track_x(), "track_x", 
				track_dim, "x particle track", "(?)");
			save_vector_1d(out_group, imp_stats.get_track_y(), "track_y", 
				track_dim, "y particle track", "(?)");
			save_vector_1d(out_group, imp_stats.get_track_z(), "track_z", 
				track_dim, "z particle track", "(?)");
			save_vector_1d(out_group, imp_stats.get_track_vx(), "track_vx", 
				track_dim, "x particle track velocity", "(?)");
			save_vector_1d(out_group, imp_stats.get_track_vy(), "track_vy", 
				track_dim, "y particle track velocity", "(?)");
			save_vector_1d(out_group, imp_stats.get_track_vz(), "track_vz", 
				track_dim, "z particle track velocity", "(?)");
		}
	}

	/**
	* Returns the netCDF type matching T: double, float or int.
	*
	* @tparam T C++ type of the values.
	* @return netCDF type.
	*/
	template <typename T>
	netCDF::NcType nc_type()
	{
		if constexpr (std::is_same_v<T, double>) return netCDF::ncDouble;
		else if constexpr (std::is_same_v<T, float>) return netCDF::ncFloat;
		else
		{
			static_assert(std::is_same_v<T, int>, "nc_type only handles "
				"double, float and int.");
			return netCDF::ncInt;
		}
	}

	/**
	* Saves one zone of a ZonedVector as a variable in nc_group.
	*
	* Only the zone's interior cells are saved, without ghost cells, as
	* (dim1, nx, ny, nz) flattened in C order. dims must match that shape,
	* without dim1 if it's 1 (e.g., geometry). An empty vector (data the
	* background doesn't have) is skipped without a warning.
	*
	* @tparam T Type of the values: double, float or int.
	* @param nc_group Group to save into.
	* @param vec Values to save.
	* @param zone Zone to save.
	* @param var_name Name of the variable.
	* @param dims Dimensions of the variable.
	* @param description Description attribute.
	* @param units Units attribute.
	*/
	template <typename T>
	void save_zone_var(const netCDF::NcGroup& nc_group,
		const Vectors::ZonedVector<T>& vec, const int zone,
		const std::string& var_name, const std::vector<netCDF::NcDim>& dims,
		const std::string& description, const std::string& units)
	{
		if (vec.get_data().empty()) return;

		std::vector<T> data {};
		data.reserve(static_cast<std::size_t>(vec.get_dim1())
			* vec.get_nx(zone) * vec.get_ny(zone) * vec.get_nz(zone));
		for (int t {}; t < vec.get_dim1(); ++t)
		for (int i {}; i < vec.get_nx(zone); ++i)
		for (int j {}; j < vec.get_ny(zone); ++j)
		for (int k {}; k < vec.get_nz(zone); ++k)
		{
			data.push_back(vec(zone, t, i, j, k));
		}

		netCDF::NcVar var {nc_group.addVar(var_name, nc_type<T>(), dims)};
		var.putVar(data.data());
		var.putAtt("description", description);
		var.putAtt("units", units);
	}

	/// A background or statistics array to save, with its name, description
	/// and units in the file
	struct ZoneVar
	{
		const Vectors::ZonedVector<BkgFPType>* vec;
		std::string name;
		std::string description;
		std::string units;
	};

	/**
	* Saves the zones of a background with more than one zone.
	*
	* A zones group holds nx, ny and nz (the number of cells in each zone)
	* and faces, what is across each face of each zone (rows xmin, xmax,
	* ymin, ymax, zmin, zmax): a zone number from 1, matching the zone
	* groups (the zone itself means periodic), -1 for an absorbing face, -2
	* for a core face or -4 for a wall face.
	*
	* Then each zone gets a group zone1, zone2, ... with its own x, y, z and
	* grid x, grid y and grid z dimensions, and geometry, background and
	* output subgroups with the same variables as a single-zone file, for
	* the zone's interior cells. The geometry also has the wall mask chi, if
	* the background has one, and the cell corners in the first toroidal
	* plane as R_corners and Z_corners (m), for plotting. Data the background
	* doesn't have (e.g., the metric coefficients of a SOLEDGE3X background)
	* is left out.
	*
	* @param nc_file File to save into. The time dimension must already
	*        exist.
	* @param opts Simulation options, used for calc_grad_elec.
	* @param bkg Background plasma to save.
	* @param imp_stats Impurity statistics to save.
	*/
	void save_zones(const netCDF::NcFile& nc_file,
		const Options::Options& opts, const Background::Background& bkg,
		Impurity::Statistics& imp_stats)
	{
		const Zones::Zones& zones {bkg.get_zones()};
		const int nzones {zones.get_nzones()};
		const netCDF::NcDim dim_time {nc_file.getDim("time")};

		// Zone table
		netCDF::NcGroup zones_group {nc_file.addGroup("zones")};
		const netCDF::NcDim dim_zone {zones_group.addDim("zone", nzones)};
		const netCDF::NcDim dim_face {zones_group.addDim("face",
			Zones::num_faces)};
		save_vector_1d(zones_group, zones.get_nx_all(), "nx", dim_zone,
			"number of x cells in each zone", "");
		save_vector_1d(zones_group, zones.get_ny_all(), "ny", dim_zone,
			"number of y cells in each zone", "");
		save_vector_1d(zones_group, zones.get_nz_all(), "nz", dim_zone,
			"number of z cells in each zone", "");
		std::vector<int> faces {};
		for (int f {}; f < Zones::num_faces; ++f)
		for (int zone {}; zone < nzones; ++zone)
		{
			const int code {zones.get_face(zone,
				static_cast<Zones::Face>(f))};
			faces.push_back(code >= 0 ? code + 1 : code);
		}
		netCDF::NcVar faces_var {zones_group.addVar("faces", netCDF::ncInt,
			{dim_face, dim_zone})};
		faces_var.putVar(faces.data());
		faces_var.putAtt("description", "what is across each face (rows "
			"xmin, xmax, ymin, ymax, zmin, zmax) of each zone: a zone number "
			"from 1 (the zone itself is periodic), -1 absorbing, -2 core or -4 "
			"wall");

		// Geometry arrays at the cell centers
		const std::vector<ZoneVar> geometry {
			{&bkg.get_X(), "X", "X cell centers", "(m)"},
			{&bkg.get_Y(), "Y", "Y cell centers", "(m)"},
			{&bkg.get_Z(), "Z", "Z cell centers", "(m)"},
			{&bkg.get_J(), "J", "Jacobian", "()"},
			{&bkg.get_gij_00(), "gij_00", "metric coefficient (00)", "(?)"},
			{&bkg.get_gij_01(), "gij_01", "metric coefficient (01)", "(?)"},
			{&bkg.get_gij_02(), "gij_02", "metric coefficient (02)", "(?)"},
			{&bkg.get_gij_11(), "gij_11", "metric coefficient (11)", "(?)"},
			{&bkg.get_gij_12(), "gij_12", "metric coefficient (12)", "(?)"},
			{&bkg.get_gij_22(), "gij_22", "metric coefficient (22)", "(?)"},
			{&bkg.get_dXdx(), "dXdx", "tangent basis vector (dX/dx)", "(?)"},
			{&bkg.get_dYdx(), "dYdx", "tangent basis vector (dY/dx)", "(?)"},
			{&bkg.get_dZdx(), "dZdx", "tangent basis vector (dZ/dx)", "(?)"},
			{&bkg.get_dXdy(), "dXdy", "tangent basis vector (dX/dy)", "(?)"},
			{&bkg.get_dYdy(), "dYdy", "tangent basis vector (dY/dy)", "(?)"},
			{&bkg.get_dZdy(), "dZdy", "tangent basis vector (dZ/dy)", "(?)"},
			{&bkg.get_dXdz(), "dXdz", "tangent basis vector (dX/dz)", "(?)"},
			{&bkg.get_dYdz(), "dYdz", "tangent basis vector (dY/dz)", "(?)"},
			{&bkg.get_dZdz(), "dZdz", "tangent basis vector (dZ/dz)", "(?)"},
			{&bkg.get_dxdX(), "dxdX", "reciprocal basis vector (dx/dX)",
				"(?)"},
			{&bkg.get_dxdY(), "dxdY", "reciprocal basis vector (dx/dY)",
				"(?)"},
			{&bkg.get_dxdZ(), "dxdZ", "reciprocal basis vector (dx/dZ)",
				"(?)"},
			{&bkg.get_dydX(), "dydX", "reciprocal basis vector (dy/dX)",
				"(?)"},
			{&bkg.get_dydY(), "dydY", "reciprocal basis vector (dy/dY)",
				"(?)"},
			{&bkg.get_dydZ(), "dydZ", "reciprocal basis vector (dy/dZ)",
				"(?)"},
			{&bkg.get_dzdX(), "dzdX", "reciprocal basis vector (dz/dX)",
				"(?)"},
			{&bkg.get_dzdY(), "dzdY", "reciprocal basis vector (dz/dY)",
				"(?)"},
			{&bkg.get_dzdZ(), "dzdZ", "reciprocal basis vector (dz/dZ)",
				"(?)"},
			{&bkg.get_b_x(), "b_x", "dual basis magnetic vector (x)", "(?)"},
			{&bkg.get_b_y(), "b_y", "dual basis magnetic vector (y)", "(?)"},
			{&bkg.get_b_z(), "b_z", "dual basis magnetic vector (z)", "(?)"},
			{&bkg.get_chi(), "chi", "wall mask (1 inside the wall)", "()"}};

		std::vector<ZoneVar> background {
			{&bkg.get_ne(), "ne", "electron density", "(m-3)"},
			{&bkg.get_te(), "Te", "electron temperature", "(eV)"},
			{&bkg.get_ti(), "Ti", "ion temperature", "(eV)"},
			{&bkg.get_vp(), "Vp", "plasma potential", "(V)"},
			{&bkg.get_eX(), "E_X", "electric field (X)", "(V/m)"},
			{&bkg.get_eY(), "E_Y", "electric field (Y)", "(V/m)"},
			{&bkg.get_eZ(), "E_Z", "electric field (Z)", "(V/m)"},
			{&bkg.get_uX(), "Ui_X", "background flow (X)", "(m/s)"},
			{&bkg.get_uY(), "Ui_Y", "background flow (Y)", "(m/s)"},
			{&bkg.get_uZ(), "Ui_Z", "background flow (Z)", "(m/s)"},
			{&bkg.get_bX(), "B_X", "magnetic field (X)", "(T)"},
			{&bkg.get_bY(), "B_Y", "magnetic field (Y)", "(T)"},
			{&bkg.get_bZ(), "B_Z", "magnetic field (Z)", "(T)"},
			{&bkg.get_bR(), "B_R", "magnetic field (R)", "(T)"},
			{&bkg.get_gradbX(), "dBdX", "magnetic field gradient (X)",
				"(T/m)"},
			{&bkg.get_gradbY(), "dBdY", "magnetic field gradient (Y)",
				"(T/m)"},
			{&bkg.get_gradbZ(), "dBdZ", "magnetic field gradient (Z)",
				"(T/m)"}};
		if (opts.calc_grad_elec_int() == 1)
		{
			background.push_back({&bkg.get_gradeX(), "dEdX",
				"electric field gradient (X)", "(V/m2)"});
			background.push_back({&bkg.get_gradeY(), "dEdY",
				"electric field gradient (Y)", "(V/m2)"});
			background.push_back({&bkg.get_gradeZ(), "dEdZ",
				"electric field gradient (Z)", "(V/m2)"});
		}

		const std::vector<ZoneVar> output {
			{&imp_stats.get_density(), "nz", "impurity density", "(m-3)"},
			{&imp_stats.get_vX(), "v_X", "average impurity X velocity",
				"(m/s)"},
			{&imp_stats.get_vY(), "v_Y", "average impurity Y velocity",
				"(m/s)"},
			{&imp_stats.get_vZ(), "v_Z", "average impurity Z velocity",
				"(m/s)"},
			{&imp_stats.get_vx(), "v_x", "average impurity x velocity",
				"(?/s)"},
			{&imp_stats.get_vy(), "v_y", "average impurity y velocity",
				"(?/s)"},
			{&imp_stats.get_vz(), "v_z", "average impurity z velocity",
				"(?/s)"},
			{&imp_stats.get_charge(), "qz", "average impurity charge", "()"},
			{&imp_stats.get_s(), "nanbu_s",
				"average Nanbu imp-ion collisionality strength", "()"}};

		for (int zone {}; zone < nzones; ++zone)
		{
			netCDF::NcGroup zone_group {nc_file.addGroup("zone"
				+ std::to_string(zone + 1))};
			const netCDF::NcDim dx {zone_group.addDim("x",
				zones.get_nx(zone))};
			const netCDF::NcDim dy {zone_group.addDim("y",
				zones.get_ny(zone))};
			const netCDF::NcDim dz {zone_group.addDim("z",
				zones.get_nz(zone))};
			const netCDF::NcDim gdx {zone_group.addDim("grid x",
				zones.get_nx(zone) + 1)};
			const netCDF::NcDim gdy {zone_group.addDim("grid y",
				zones.get_ny(zone) + 1)};
			const netCDF::NcDim gdz {zone_group.addDim("grid z",
				zones.get_nz(zone) + 1)};

			// Geometry
			netCDF::NcGroup geo {zone_group.addGroup("geometry")};
			save_vector_1d(geo, zones.get_x(zone), "x", dx,
				"x cell centers", "(?)");
			save_vector_1d(geo, zones.get_y(zone), "y", dy,
				"y cell centers", "(?)");
			save_vector_1d(geo, zones.get_z(zone), "z", dz,
				"z cell centers", "(?)");
			save_vector_1d(geo, zones.get_grid_x(zone), "grid_x", gdx,
				"x grid edges", "(?)");
			save_vector_1d(geo, zones.get_grid_y(zone), "grid_y", gdy,
				"y grid edges", "(?)");
			save_vector_1d(geo, zones.get_grid_z(zone), "grid_z", gdz,
				"z grid edges", "(?)");
			for (const auto& v : geometry)
			{
				save_zone_var(geo, *v.vec, zone, v.name, {dx, dy, dz},
					v.description, v.units);
			}
			save_zone_var(geo, bkg.get_grid_X(), zone, "grid_X",
				{gdx, gdy, gdz}, "X grid edges", "(m)");
			save_zone_var(geo, bkg.get_grid_Y(), zone, "grid_Y",
				{gdx, gdy, gdz}, "Y grid edges", "(m)");
			save_zone_var(geo, bkg.get_grid_Z(), zone, "grid_Z",
				{gdx, gdy, gdz}, "Z grid edges", "(m)");

			// Poloidal cell corners in the first toroidal plane
			const auto& grid_X {bkg.get_grid_X()};
			if (!grid_X.get_data().empty())
			{
				std::vector<double> R {}, Z {};
				for (int i {}; i <= zones.get_nx(zone); ++i)
				for (int j {}; j <= zones.get_ny(zone); ++j)
				{
					R.push_back(std::hypot(grid_X(zone, i, j, 0),
						bkg.get_grid_Y()(zone, i, j, 0)));
					Z.push_back(bkg.get_grid_Z()(zone, i, j, 0));
				}
				netCDF::NcVar R_var {geo.addVar("R_corners",
					netCDF::ncDouble, {gdx, gdy})};
				R_var.putVar(R.data());
				R_var.putAtt("description", "R of the cell corners");
				R_var.putAtt("units", "(m)");
				netCDF::NcVar Z_var {geo.addVar("Z_corners",
					netCDF::ncDouble, {gdx, gdy})};
				Z_var.putVar(Z.data());
				Z_var.putAtt("description", "Z of the cell corners");
				Z_var.putAtt("units", "(m)");
			}

			// Background
			netCDF::NcGroup bkg_group {zone_group.addGroup("background")};
			for (const auto& v : background)
			{
				save_zone_var(bkg_group, *v.vec, zone, v.name,
					{dim_time, dx, dy, dz}, v.description, v.units);
			}

			// Output
			netCDF::NcGroup out_group {zone_group.addGroup("output")};
			save_zone_var(out_group, imp_stats.get_counts(), zone, "Nz",
				{dim_time, dx, dy, dz}, "impurity counts", "");
			for (const auto& v : output)
			{
				save_zone_var(out_group, *v.vec, zone, v.name,
					{dim_time, dx, dy, dz}, v.description, v.units);
			}
		}
	}

	/**
	* Saves everything into a netCDF file named <case_name>.nc.
	*
	* An existing file with that name is deleted first. The dimensions are
	* created here, then the commit hash, input options, geometry,
	* background and impurity results are saved (see the @file block for
	* the layout).
	*
	* @param bkg Background plasma, used for the dimensions and saved in
	*        the geometry and background groups.
	* @param imp_stats Impurity statistics to save.
	* @param opts Simulation options.
	*/
	void save_netcdf(const Background::Background& bkg, 
		Impurity::Statistics& imp_stats, const Options::Options& opts)
	{
		std::cout << "Saving netCDF file...\n";

		// Reusable string to pass in descriptions
		std::string desc {};

		// Open new netCDF file named after this case. replace tells netCDF
		// to overwrite the file if it exists. NcFile can't be copied, so it
		// is passed to the save functions by reference, and it closes
		// automatically once it goes out of scope.
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
		// First dimension is time, then x, y, z. With more than one zone,
		// each zone's group has its own x, y and z (see save_zones).
		netCDF::NcDim dim1 {nc_file.addDim("time", bkg.get_dim1())};
		const bool one_zone {bkg.get_nzones() == 1};
		if (one_zone)
		{
			nc_file.addDim("x", bkg.get_dim2());
			nc_file.addDim("y", bkg.get_dim3());
			nc_file.addDim("z", bkg.get_dim4());

			// Dimensions for the grid points. They're just bigger than the
			// cell center lengths.
			nc_file.addDim("grid x", bkg.get_dim2() + 1);
			nc_file.addDim("grid y", bkg.get_dim3() + 1);
			nc_file.addDim("grid z", bkg.get_dim4() + 1);

			// Dimensions for the particle track
			nc_file.addDim("track", imp_stats.get_track_x().size());
		}

		// Details about the simulation (system, time run, time spent, etc.)
		// To-do

		// Save commit hash.
		const std::string commit_hash {GIT_COMMIT_HASH};
		desc = "commit hash";
		save_string(nc_file, commit_hash, "commit", dim_str, desc);

		// Save all input options used for this case
		save_input_options(nc_file, opts);

		if (one_zone)
		{
			// Save all geometry related variables
			save_geometry(nc_file, bkg);

			// Save all background related variables
			save_background(nc_file, opts, bkg);

			// Save impurity transport results
			save_output(nc_file, opts, bkg, imp_stats);
		}

		// With more than one zone, only the times are in the top-level
		// geometry group, and everything else is per zone
		else
		{
			netCDF::NcGroup geo_group {nc_file.addGroup("geometry")};
			save_vector_1d(geo_group, bkg.get_times(), "time", dim1,
				"time for each frame", "(s)");
			save_zones(nc_file, opts, bkg, imp_stats);
		}

		std::cout << "Saved as " << nc_filename << '\n';
	}

}  // namespace SaveResults
