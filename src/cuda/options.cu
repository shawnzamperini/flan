#include <cuda_runtime.h>

#include "options.h"
#include "options_device.h"


namespace Options
{
	// Copy select options to device, returning a OptionsDevice object with
	// pointers to memory location on device.
	OptionsDevice* Options::to_device(int device_id)
	{
		cudaSetDevice(device_id);
		
		// Create a OptionsDevice on the host, copy all the options we need
		// into it, and then copy the struct to the device at the end and
		// return the struct that can be copied into kernels
		OptionsDevice opts_h {};

		opts_h.tstart_opt_int = m_imp_tstart_opt_int;
		opts_h.tstart_val     = m_imp_tstart_val;
		opts_h.trange_min     = m_imp_trange_min;
		opts_h.trange_max     = m_imp_trange_max;

		opts_h.xstart_opt_int = m_imp_xstart_opt_int;
		opts_h.xstart_val     = m_imp_xstart_val;
		opts_h.xrange_min     = m_imp_xrange_min;
		opts_h.xrange_max     = m_imp_xrange_max;

		opts_h.ystart_opt_int = m_imp_ystart_opt_int;
		opts_h.ystart_val     = m_imp_ystart_val;
		opts_h.yrange_min     = m_imp_yrange_min;
		opts_h.yrange_max     = m_imp_yrange_max;

		opts_h.zstart_opt_int = m_imp_zstart_opt_int;
		opts_h.zstart_val     = m_imp_zstart_val;
		opts_h.zrange_min     = m_imp_zrange_min;
		opts_h.zrange_max     = m_imp_zrange_max;

		opts_h.init_charge    = m_imp_init_charge;

		OptionsDevice* opts_d = nullptr;
		cudaMalloc(&opts_d, sizeof(OptionsDevice));
		cudaMemcpy(opts_d, &opts_h, sizeof(OptionsDevice), cudaMemcpyHostToDevice);

		return opts_d;
	}

	// Free the device-side struct
	void free_opts(OptionsDevice* opts_d, int device_id)
	{
		cudaSetDevice(device_id);
		cudaFree(opts_d);
	}


}  // namespace Options
