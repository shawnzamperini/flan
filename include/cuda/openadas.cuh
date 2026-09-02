#include "background_device.h"
#include "openadas_device.h"
#include "pcg32.h"
#include "slots_device.h"


namespace OpenADAS
{
	/**
	* @brief Wrapper to call ioniz_recomb_kernel, which updates particle
	*   charge based on ionization/recombination rates.
	*/
	void ioniz_recomb_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d,
		const OpenADASDevice& oa_ioniz_d, 
		const OpenADASDevice& oa_recomb_d, const double dt, 
		int& ioniz_warnings, int& recomb_warnings, pcg32* rngs_d);

}  // OpenADAS
