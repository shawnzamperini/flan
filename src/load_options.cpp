#include <iostream>

#include "flan_types.h"  // For Inputs type
#include "load_options.h"
#include "options.h"

namespace Options
{
	Options load_options(const Inputs& inpts)
	{
		
		// First load an Options with all the defaults preloaded
		//Options opts {load_defaults()};

		// Create Options with default values
		Options opts {};

		// Then load in the user input, overwriting with any options they
		// provided
		load_input(opts, inpts);

		return opts;
	}

	void load_input(Options& opts, const Inputs& inpts)
	{
		std::cout << "Loading input...\n";
		
		// Loop through every input option. key is a string of the option, and
		// var stands for "variant", which is a type that can contain multiple
		// types. We dissect var into the underlying type within assign_option.
		for (const auto& [key, var] : inpts)
		{
			// Overwrite default with given value. Types must match. This is a
			// little tricky, but it's done this way to handle multiple types.
			// First we define set_func, which binds opts to matching
			// assignment function from the Options class. This is required
			// because C++ does not allowing passing member functions, you
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
				auto set_func = std::bind(&Options::set_gkyl_elec_mass_amu, &opts, 
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}
			
			else if (key == "gkyl_ion_mass_amu") 
			{
				auto set_func = std::bind(&Options::set_gkyl_ion_mass_amu, &opts, 
					std::placeholders::_1);
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

			else if (key == "min_xbound_type") 
			{
				auto set_func = std::bind(&Options::set_min_xbound_type, &opts, 
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}
			
			else if (key == "lcfs_x") 
			{
				auto set_func = std::bind(&Options::set_lcfs_x, &opts, 
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
				auto set_func = std::bind(&Options::set_imp_time_step_min, &opts, 
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}
			
			else if (key == "imp_source_scale_fact") 
			{
				auto set_func = std::bind(&Options::set_imp_source_scale_fact, 
					&opts, std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}
			
			else if (key == "imp_vel_stats") 
			{
				auto set_func = std::bind(&Options::set_imp_vel_stats, &opts, 
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
			}
			
			else if (key == "imp_xbound_buffer") 
			{
				auto set_func = std::bind(&Options::set_imp_xbound_buffer, &opts, 
					std::placeholders::_1);
				assign_option<double>(set_func, var, key);
			}
			
			else if (key == "imp_iz_recomb") 
			{
				auto set_func = std::bind(&Options::set_imp_iz_recomb, &opts, 
					std::placeholders::_1);
				assign_option<std::string>(set_func, var, key);
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
			/* Hand-written one, don't want to delete until I check the above.	
			if (key == "imp_num") 
			{
				auto set_func = std::bind(&Options::set_imp_num, &opts, 
					std::placeholders::_1);
				assign_option<int>(set_func, var, key);
			}
			*/
			else if (key == "mapc2p") 
			{
				auto set_func = std::bind(&Options::set_mapc2p, &opts, 
					std::placeholders::_1);
				assign_option<Mapc2p_ptr>(set_func, var, key);
			}
			else
			{
				std::cerr << "Error! " << key << " is not a recognized "
					<< "option\n";
			}
		}
	}

	// Returns true if the variant (var) is of type T, false if not.
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
			// Print error message, letting user know what type the variable
			// should be.
			std::cerr << "Error! " << var_name << " is the wrong type. "
				<< "Expected ";
			if (std::is_same<T, int>::value) std::cerr << "int";
			else if (std::is_same<T, double>::value) std::cerr << "double";
			else if (std::is_same<T, std::string>::value) std::cerr 
				<< "string";
			else if (std::is_same<T, Mapc2p_ptr>::value) std::cerr 
				<< "Mapc2p_ptr";
			std::cerr << ".\n";

			return false;
		}
	}

	// Uses a set function (set_func) from our Options object to 
	// set the value stored in the variant (var). If the data in var is not of
	// type T, then the value is not stored and an error is printed. 
	template <typename T>
	void assign_option(std::function<void(T)> set_func, auto var, 
		const std::string_view var_name)
	{
		// Call set function only if the data in var is of the correct type T.
		if (var_correct_type<T>(var, var_name))
		{
			set_func(std::get<T>(var));
		}
		else
		{
			std::cerr << var_name << " not set.\n";
		}
	}

	// Check that various options are valid with each other
	void crosscheck_options(Options& opts)
	{
		// Particle splitting based on ionization/recombination requires
		// ioniz_recomb to be true

		// Variable time step is based on collisions, so requires it to be on

		// Nanbu collisions require BiMaxwellian moments from Gkeyll
		if (opts.imp_collisions() == "nanbu" 
			&& opts.gkyl_moment_type() != "bimaxwellian")
		{
			std::cout << "Warning! Nanbu collisions requires BiMaxwellianMoment"
				<< " files from Gkeyll. Collisions are being turned OFF.\n";
			opts.set_imp_collisions("off");
		}
	}

}
