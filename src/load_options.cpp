/**
* @file  load_options.cpp
*
* @brief Routines to handle loading input options into an Options object.
*
* The input file assembles the options into an Inputs object, which maps
* each option name to its value. Each value is checked to be the right type
* and then passed to the matching setter in Options. Options not given in
* the input file keep their default values. See the Input Options page of the
* documentation for every option.
*/

#include <iostream>
#include <sstream>
#include <string>

#include "errors.h"
#include "flan_types.h"  // For Inputs type
#include "load_options.h"
#include "options.h"


namespace Options
{
	/**
	* Creates an Options object from the input options.
	*
	* Starts from an Options object with the default values, then
	* overwrites them with any options given in inpts. Finally the options
	* are checked to be compatible with each other with crosscheck_options,
	* which may change some of them.
	*
	* @param inpts Input options assembled in a simulation's input file.
	* @return Options object with the values from inpts. Options that are
	*         not in inpts keep their default values.
	*/
	Options load_options(const Inputs& inpts)
	{
		// Create Options with default values
		Options opts {};

		// Then load in the user input, overwriting with any options they
		// provided
		load_input(opts, inpts);

		// Check that the options are compatible with each other, changing
		// any that aren't
		crosscheck_options(opts);

		return opts;
	}

	/**
	* Loads each input option into an Options object, overwriting the
	* default value.
	*
	* Each option name is matched to its setter in Options and the value is
	* set with assign_option. The run is stopped with Errors::fatal_error
	* for an option name that is not recognized (e.g., a typo) or a value
	* of the wrong type, since running with the default instead would
	* silently ignore the option.
	*
	* @param opts Options object to load the options into. Generally should
	*        be passed in right after constructing it.
	* @param inpts Input options assembled in a simulation's input file.
	*/
	void load_input(Options& opts, const Inputs& inpts)
	{
		// Loop through every input option. key is a string of the option, and
		// var stands for "variant", which is a type that can contain multiple
		// types. We dissect var into the underlying type within assign_option.
		for (const auto& [key, var] : inpts)
		{
			// Overwrite default with given value. Types must match. This is a
			// little tricky, but it's done this way to handle multiple types.
			// First we define set_func, which binds opts to matching
			// assignment function from the Options class. This is required
			// because C++ does not allow passing member functions, you
			// need to bind them first. placeholders::_1 is there to say the
			// function still needs to be passed in an argument to be used -
			// this will be the value to set.
			// Then we pass the set_func along with the variant (var) and
			// name of the variable (key) to assign_option, making sure to
			// use the correct template type for the argument. The value is
			// set within assign_option.
			if (key == "case_name")
			{
				auto set_func = std::bind(&Options::set_case_name, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "use_gpu")
			{
				auto set_func = std::bind(&Options::set_use_gpu, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "seed")
			{
				auto set_func = std::bind(&Options::set_seed, &opts,
					std::placeholders::_1);
				assign_option<int>(set_func, var, key);
			}

			else if (key == "slot_cap")
			{
				auto set_func = std::bind(&Options::set_slot_cap, &opts,
					std::placeholders::_1);
				assign_option<int>(set_func, var, key);
			}

			else if (key == "bkg_source")
			{
				auto set_func = std::bind(&Options::set_bkg_source, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "test_opt")
			{
				auto set_func = std::bind(&Options::set_test_opt, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "save_track")
			{
				auto set_func = std::bind(&Options::set_save_track, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "gkyl_dir")
			{
				auto set_func = std::bind(&Options::set_gkyl_dir, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "gkyl_casename")
			{
				auto set_func = std::bind(&Options::set_gkyl_casename, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "gkyl_frame_start")
			{
				auto set_func = std::bind(&Options::set_gkyl_frame_start, &opts,
					std::placeholders::_1);
				assign_option<int>(set_func, var, key);
			}

			else if (key == "gkyl_frame_end")
			{
				auto set_func = std::bind(&Options::set_gkyl_frame_end, &opts,
					std::placeholders::_1);
				assign_option<int>(set_func, var, key);
			}

			else if (key == "gkyl_elec_name")
			{
				auto set_func = std::bind(&Options::set_gkyl_elec_name, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "gkyl_ion_name")
			{
				auto set_func = std::bind(&Options::set_gkyl_ion_name, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "gkyl_elec_mass_amu")
			{
				auto set_func = std::bind(&Options::set_gkyl_elec_mass_amu,
					&opts, std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "gkyl_ion_mass_amu")
			{
				auto set_func = std::bind(&Options::set_gkyl_ion_mass_amu,
					&opts, std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "gkyl_file_type")
			{
				auto set_func = std::bind(&Options::set_gkyl_file_type, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "gkyl_moment_type")
			{
				auto set_func = std::bind(&Options::set_gkyl_moment_type, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "calc_grad_elec")
			{
				auto set_func = std::bind(&Options::set_calc_grad_elec, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "min_ne")
			{
				auto set_func = std::bind(&Options::set_min_ne, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "min_te")
			{
				auto set_func = std::bind(&Options::set_min_te, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "min_ti")
			{
				auto set_func = std::bind(&Options::set_min_ti, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "tbound_type")
			{
				auto set_func = std::bind(&Options::set_tbound_type, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "min_xbound_type")
			{
				auto set_func = std::bind(&Options::set_min_xbound_type, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "max_xbound_type")
			{
				auto set_func = std::bind(&Options::set_max_xbound_type, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "min_ybound_type")
			{
				auto set_func = std::bind(&Options::set_min_ybound_type, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "max_ybound_type")
			{
				auto set_func = std::bind(&Options::set_max_ybound_type, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "min_zbound_type")
			{
				auto set_func = std::bind(&Options::set_min_zbound_type, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "max_zbound_type")
			{
				auto set_func = std::bind(&Options::set_max_zbound_type, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "lcfs_x")
			{
				auto set_func = std::bind(&Options::set_lcfs_x, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "lcfs_y")
			{
				auto set_func = std::bind(&Options::set_lcfs_y, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "lcfs_z")
			{
				auto set_func = std::bind(&Options::set_lcfs_z, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "sep_x_bc_xp_z1")
			{
				auto set_func = std::bind(&Options::set_sep_x_bc_xp_z1, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "sep_x_bc_xp_z2")
			{
				auto set_func = std::bind(&Options::set_sep_x_bc_xp_z2, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_atom_num")
			{
				auto set_func = std::bind(&Options::set_imp_atom_num, &opts,
					std::placeholders::_1);
				assign_option<int>(set_func, var, key);
			}

			else if (key == "imp_mass_amu")
			{
				auto set_func = std::bind(&Options::set_imp_mass_amu, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_init_charge")
			{
				auto set_func = std::bind(&Options::set_imp_init_charge, &opts,
					std::placeholders::_1);
				assign_option<int>(set_func, var, key);
			}

			else if (key == "imp_num")
			{
				auto set_func = std::bind(&Options::set_imp_num, &opts,
					std::placeholders::_1);
				assign_option<int>(set_func, var, key);
			}

			else if (key == "imp_tstart_opt")
			{
				auto set_func = std::bind(&Options::set_imp_tstart_opt, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "imp_tstart_val")
			{
				auto set_func = std::bind(&Options::set_imp_tstart_val, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_trange_min")
			{
				auto set_func = std::bind(&Options::set_imp_trange_min, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_trange_max")
			{
				auto set_func = std::bind(&Options::set_imp_trange_max, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_xstart_opt")
			{
				auto set_func = std::bind(&Options::set_imp_xstart_opt, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "imp_xstart_val")
			{
				auto set_func = std::bind(&Options::set_imp_xstart_val, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_xrange_min")
			{
				auto set_func = std::bind(&Options::set_imp_xrange_min, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_xrange_max")
			{
				auto set_func = std::bind(&Options::set_imp_xrange_max, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_ystart_opt")
			{
				auto set_func = std::bind(&Options::set_imp_ystart_opt, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "imp_ystart_val")
			{
				auto set_func = std::bind(&Options::set_imp_ystart_val, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_yrange_min")
			{
				auto set_func = std::bind(&Options::set_imp_yrange_min, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_yrange_max")
			{
				auto set_func = std::bind(&Options::set_imp_yrange_max, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_zstart_opt")
			{
				auto set_func = std::bind(&Options::set_imp_zstart_opt, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "imp_zstart_val")
			{
				auto set_func = std::bind(&Options::set_imp_zstart_val, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_zrange_min")
			{
				auto set_func = std::bind(&Options::set_imp_zrange_min, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_zrange_max")
			{
				auto set_func = std::bind(&Options::set_imp_zrange_max, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_temp_start_opt")
			{
				auto set_func = std::bind(&Options::set_imp_temp_start_opt,
					&opts, std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "imp_temp_start_val")
			{
				auto set_func = std::bind(&Options::set_imp_temp_start_val,
					&opts, std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_collisions")
			{
				auto set_func = std::bind(&Options::set_imp_collisions, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "var_red_split")
			{
				auto set_func = std::bind(&Options::set_var_red_split, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "var_red_import")
			{
				auto set_func = std::bind(&Options::set_var_red_import, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "var_red_freq")
			{
				auto set_func = std::bind(&Options::set_var_red_freq,
					&opts, std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "var_red_min_weight")
			{
				auto set_func = std::bind(
					&Options::set_var_red_min_weight, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "var_red_med_mod")
			{
				auto set_func = std::bind(
					&Options::set_var_red_med_mod, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "var_red_rusrol")
			{
				auto set_func = std::bind(&Options::set_var_red_rusrol, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "var_red_rusrol_prob")
			{
				auto set_func = std::bind(&Options::set_var_red_rusrol_prob,
					&opts, std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_time_step_opt")
			{
				auto set_func = std::bind(&Options::set_imp_time_step_opt,
					&opts, std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "imp_time_step")
			{
				auto set_func = std::bind(&Options::set_imp_time_step, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_time_step_min")
			{
				auto set_func = std::bind(&Options::set_imp_time_step_min,
					&opts, std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_source_scale_fact")
			{
				auto set_func = std::bind(&Options::set_imp_source_scale_fact,
					&opts, std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_xbound_buffer")
			{
				auto set_func = std::bind(&Options::set_imp_xbound_buffer,
					&opts, std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_ybound_buffer")
			{
				auto set_func = std::bind(&Options::set_imp_ybound_buffer,
					&opts, std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_zbound_buffer")
			{
				auto set_func = std::bind(&Options::set_imp_zbound_buffer,
					&opts, std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "imp_iz_recomb")
			{
				auto set_func = std::bind(&Options::set_imp_iz_recomb, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "imp_max_time")
			{
				auto set_func = std::bind(&Options::set_imp_max_time, &opts,
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}

			else if (key == "print_interval")
			{
				auto set_func = std::bind(&Options::set_print_interval, &opts,
					std::placeholders::_1);
				assign_option<int>(set_func, var, key);
			}

			else if (key == "openadas_root")
			{
				auto set_func = std::bind(&Options::set_openadas_root, &opts,
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}

			else if (key == "openadas_year")
			{
				auto set_func = std::bind(&Options::set_openadas_year, &opts,
					std::placeholders::_1);
				assign_option<int>(set_func, var, key);
			}

			else
			{
				Errors::fatal_error(std::string(key) 
					+ " is not a recognized option.");
			}
		}
	}

	/**
	* Checks if the value in a variant is of the expected type.
	*
	* If not, the run is stopped with Errors::fatal_error, saying which type
	* was expected.
	*
	* @tparam T Expected type of the value: int, double, std::string or
	*         Mapc2p_ptr.
	* @param var Variant holding the value of an input option. See Inputs in
	*        flan_types.h.
	* @param var_name Name of the input option, used in the error message.
	* @return True if the value in var is of type T. Doesn't return if
	*         not.
	*/
	template <typename T>
	bool var_correct_type(auto var, const std::string_view var_name)
	{
		// get_if<T>(&var) returns a pointer of the underlying type of the
		// variant (int, double, ...) if it matches T. Otherwise, it returns
		// a nullptr. So if we call this function template with an int, and
		// the underlying type is an int then it returns that int. If the
		// underlying type is anything else, an error message is output.
		if (const T* var_ptr {std::get_if<T>(&var)})
		{
			return true;
		}
		else
		{
			// Stop the run, letting user know what type the variable
			// should be.
			std::ostringstream msg {};
			msg << var_name << " is the wrong type. Expected ";
			if (std::is_same<T, int>::value) msg << "int";
			else if (std::is_same<T, double>::value) msg << "double";
			else if (std::is_same<T, std::string>::value) msg << "string";
			else if (std::is_same<T, Mapc2p_ptr>::value) msg << "Mapc2p_ptr";
			msg << ".";
			Errors::fatal_error(msg.str());
		}
	}

	/**
	* Sets an option to the value stored in a variant, overwriting the
	* default.
	*
	* If the value in var is not of type T, var_correct_type stops the run.
	*
	* @tparam T Expected type of the value: int, double, std::string or
	*         Mapc2p_ptr.
	* @param set_func Setter from the Options object, bound to the object
	*        with std::bind.
	* @param var Variant holding the value of the input option.
	* @param var_name Name of the input option, used in the error messages.
	*/
	template <typename T>
	void assign_option(std::function<void(T)> set_func, auto var,
		const std::string_view var_name)
	{
		// Call set function only if the data in var is of the correct type
		// T. If it isn't, var_correct_type stops the run.
		if (var_correct_type<T>(var, var_name))
		{
			set_func(std::get<T>(var));
		}
	}

	/**
	* Checks that the options are compatible with each other, changing them
	* if needed.
	*
	* Currently this checks that:
	*  - Nanbu collisions are used with BiMaxwellian moments from Gkeyll.
	*    If not, collisions are turned off.
	*  - A periodic x, y or z boundary is periodic on both the minimum and
	*    maximum side. If only one side is, both are set to periodic.
	*  - GPU acceleration is off in a CPU-only build.
	*
	* A warning is printed whenever an option is changed. Placeholders are
	* left for checks of particle splitting and the variable time step.
	*
	* @param opts Options object to check and modify in place.
	*/
	void crosscheck_options(Options& opts)
	{
		// Nanbu collisions require BiMaxwellian moments from Gkeyll
		if (opts.imp_collisions() == "nanbu"
			&& opts.gkyl_moment_type() != "bimaxwellian")
		{
			std::cout << "Warning! Nanbu collisions requires BiMaxwellianMoment"
				<< " files from Gkeyll. Collisions are being turned OFF.\n";
			opts.set_imp_collisions("off");
		}

		// A periodic boundary wraps the particle to the other side of the
		// domain, so it only makes sense if both sides of a direction are
		// periodic. periodic_bc wraps at both ends and the minimum side is
		// checked first, so, e.g., a periodic min x with an absorbing max x
		// would wrap particles past max x before they could be absorbed. If
		// only one side is periodic, make both sides periodic.
		auto check_periodic_pair = [&opts](const std::string& dir,
			const std::string& min_type, const std::string& max_type,
			void (Options::*set_min)(std::string),
			void (Options::*set_max)(std::string))
		{
			if ((min_type == "periodic") != (max_type == "periodic"))
			{
				std::cout << "Warning! min_" << dir << "bound_type = "
					<< min_type << " and max_" << dir << "bound_type = "
					<< max_type << ", but a periodic boundary must be "
					<< "periodic on both sides. Setting both to periodic.\n";
				(opts.*set_min)("periodic");
				(opts.*set_max)("periodic");
			}
		};
		check_periodic_pair("x", opts.min_xbound_type(),
			opts.max_xbound_type(), &Options::set_min_xbound_type,
			&Options::set_max_xbound_type);
		check_periodic_pair("y", opts.min_ybound_type(),
			opts.max_ybound_type(), &Options::set_min_ybound_type,
			&Options::set_max_ybound_type);
		check_periodic_pair("z", opts.min_zbound_type(),
			opts.max_zbound_type(), &Options::set_min_zbound_type,
			&Options::set_max_zbound_type);

#ifndef USE_CUDA
		// Turn off GPU option if it was turned on for a CPU-only build.
		if (opts.use_gpu_int() > 0)
		{
			std::cout << "GPU acceleration requested (use_gpu = " 
				<< opts.use_gpu() << ") but Flan was built without GPU "
				"support. Turning off.";
			opts.set_use_gpu("off");
		}
#endif

	}
}
