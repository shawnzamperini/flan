/**
* @file  options_device.h
* @brief Contains struct to hold the options used in kernels for
*        device-side access.
*/
#pragma once

namespace Options
{
	/**
	* Device-side copy of the options that kernels need, made by
	* Options::to_device and freed with free_opts.
	*
	* Only a subset of the options is here. Kernels take a pointer to a
	* single copy in GPU memory (opts_d). A field added here must also be
	* set in Options::to_device. The *_opt_int fields are the integer
	* versions of string options, set by the matching Options setters.
	*/
	struct OptionsDevice
	{
		int test_opt_int;		///< test_opt: 0 gyrate, 1 exb, 2 gradb,
								///< 3 polarization, 4 curvature,
								///< 5 friction_force
		int bkg_source_int;		///< bkg_source: 0 test, 1 gkeyll

		double imp_mass_amu;	///< Impurity mass (amu)
		int imp_num;			///< Number of impurities to follow

		double min_ne;			///< Minimum electron density (m^-3)
		double min_te;			///< Minimum electron temperature (eV)
		double min_ti;			///< Minimum ion temperature (eV)

		// Starting time and position. Each opt_int is 0 for single_value
		// (start at *_val), 1 for range (uniform between *range_min and
		// *range_max) and 2 for full_range (uniform over the whole
		// simulation).
		int tstart_opt_int;		///< imp_tstart_opt
		double tstart_val;		///< imp_tstart_val (s)
		double trange_min;		///< imp_trange_min (s)
		double trange_max;		///< imp_trange_max (s)

		int xstart_opt_int;		///< imp_xstart_opt
		double xstart_val;		///< imp_xstart_val
		double xrange_min;		///< imp_xrange_min
		double xrange_max;		///< imp_xrange_max

		int ystart_opt_int;		///< imp_ystart_opt
		double ystart_val;		///< imp_ystart_val
		double yrange_min;		///< imp_yrange_min
		double yrange_max;		///< imp_yrange_max

		int zstart_opt_int;		///< imp_zstart_opt
		double zstart_val;		///< imp_zstart_val
		double zrange_min;		///< imp_zrange_min
		double zrange_max;		///< imp_zrange_max

		double imp_max_time;	///< Maximum time a particle is followed (s)

		// Copied, but not read by any kernel yet
		double lcfs_x;			///< x coordinate of the last closed flux
								///< surface
		double lcfs_y;			///< y coordinate of the last closed flux
								///< surface
		double lcfs_z;			///< z coordinate of the last closed flux
								///< surface

		int init_charge;		///< Charge state particles are launched with

		int imp_temp_start_opt_int;	///< imp_temp_start_opt: 0 single_value,
									///< 1 main_ion
		double imp_temp_start_val;	///< Starting temperature for
									///< single_value (eV)

	}; // OptionsDevice
}  // namespace Options
