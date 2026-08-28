#pragma once 

namespace Options
{
	struct OptionsDevice
	{
		// Only a subset of the options that are actually used in kernels
		int tstart_opt_int;
		double tstart_val;
		double trange_min;
		double trange_max;
		
		int xstart_opt_int;
		double xstart_val;
		double xrange_min;
		double xrange_max;

		int ystart_opt_int;
		double ystart_val;
		double yrange_min;
		double yrange_max;

		int zstart_opt_int;
		double zstart_val;
		double zrange_min;
		double zrange_max;

		int init_charge;

	}; // OptionsDevice
}  // namespace Options
