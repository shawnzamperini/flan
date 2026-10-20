/**
* @file  options.cu
*
* @brief GPU definitions of Options::to_device and free_opts.
*
* Kernels only need a few of the options, so just those are copied into an
* OptionsDevice struct in GPU memory. Kernels take a pointer to it
* (opts_d). CPU-only builds use the placeholder versions in options.cpp
* instead.
*/
#include <cuda_runtime.h>

#include "options.h"
#include "options_device.h"


namespace Options
{
	/**
	* Copies the options used in kernels to a GPU.
	*
	* The options are gathered into an OptionsDevice on the host, which is
	* then copied as a whole to a single OptionsDevice allocated on the
	* GPU. Options added to OptionsDevice must also be set here, or the
	* kernels will see zero.
	*
	* @param device_id ID of the GPU to copy to.
	* @return Device pointer to the OptionsDevice. It can be passed straight
	*         to kernels, but not dereferenced on the host. Free it with
	*         free_opts.
	*/
	OptionsDevice* Options::to_device(int device_id)
	{
		cudaSetDevice(device_id);
		
		// Create an OptionsDevice on the host, copy all the options we need
		// into it, and then copy the struct to the device at the end and
		// return the pointer to it, which kernels take
		OptionsDevice opts_h {};

		opts_h.test_opt_int   = m_test_opt_int;
		opts_h.bkg_source_int = m_bkg_source_int;

		opts_h.min_ne         = m_min_ne;
		opts_h.min_te         = m_min_te;
		opts_h.min_ti         = m_min_ti;

		opts_h.imp_mass_amu   = m_imp_mass_amu;

		opts_h.imp_num        = m_imp_num;

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

		opts_h.imp_max_time   = m_imp_max_time;

		opts_h.lcfs_x         = m_lcfs_x;
		opts_h.lcfs_y         = m_lcfs_y;
		opts_h.lcfs_z         = m_lcfs_z;

		opts_h.init_charge    = m_imp_init_charge;

		opts_h.imp_temp_start_opt_int = m_imp_temp_start_opt_int;
		opts_h.imp_temp_start_val = m_imp_temp_start_val;

		OptionsDevice* opts_d = nullptr;
		cudaMalloc(&opts_d, sizeof(OptionsDevice));
		cudaMemcpy(opts_d, &opts_h, sizeof(OptionsDevice), cudaMemcpyHostToDevice);

		return opts_d;
	}  // to_device

	/**
	* Frees the OptionsDevice allocated by Options::to_device.
	*
	* @param opts_d Device pointer returned by Options::to_device. It is
	*        not set to nullptr, since it's passed by value.
	* @param device_id ID of the GPU the options are on.
	*/
	void free_opts(OptionsDevice* opts_d, int device_id)
	{
		cudaSetDevice(device_id);
		cudaFree(opts_d);
	}  // free_opts


}  // namespace Options
