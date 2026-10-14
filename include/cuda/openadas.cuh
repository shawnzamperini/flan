#include "background_device.h"
#include "openadas_device.h"
#include "pcg32.h"
#include "slots_device.h"


namespace OpenADAS
{
	/**
	* @brief Allocates the ionization/recombination warning counters on a
	*   GPU and sets them to zero.
	*
	* Called once per GPU at the start of a run. The counters then
	* accumulate over every call to ioniz_recomb_gpu.
	*
	* @param device_id ID of the GPU to allocate the counters on.
	* @return Struct holding the device pointers to the two counters.
	*/
	IonizRecombWarningsDevice alloc_warnings(int device_id);

	/**
	* @brief Copies the warning counts from the GPU and adds them to the
	*   host-side counters.
	*
	* Blocks until all previously launched kernels have finished.
	*
	* @param warnings_d Device-side warning counters.
	* @param ioniz_warnings Host-side ionization warning counter.
	*                       Incremented by the device count.
	* @param recomb_warnings Host-side recombination warning counter.
	*                        Incremented by the device count.
	*/
	void add_warnings_to_host(const IonizRecombWarningsDevice& warnings_d,
		int& ioniz_warnings, int& recomb_warnings);

	/**
	* @brief Frees the device-side warning counters.
	*
	* @param warnings_d Device-side warning counters. The pointers are set
	*                   to nullptr.
	*/
	void free_warnings(IonizRecombWarningsDevice& warnings_d);

	/**
	* @brief Wrapper to call ioniz_recomb_kernel, which updates particle
	*   charge based on ionization/recombination rates.
	*
	* This is only a kernel launch, so it returns before the kernel has
	* finished. The number of particles whose ionization or recombination
	* probability exceeded 1 is added to the counters in warnings_d.
	*
	* @param slots_d Device-side particle slots. Charges updated in place.
	* @param bkg_d Device-side background plasma.
	* @param oa_ioniz_d Device-side OpenADAS ionization (scd) rates.
	* @param oa_recomb_d Device-side OpenADAS recombination (acd) rates.
	* @param dt Time step (s).
	* @param warnings_d Device-side warning counters from alloc_warnings.
	* @param rngs_d Device-side random number generators, one per slot.
	*/
	void ioniz_recomb_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d,
		const OpenADASDevice& oa_ioniz_d, 
		const OpenADASDevice& oa_recomb_d, const double dt, 
		const IonizRecombWarningsDevice& warnings_d, pcg32* rngs_d);

}  // OpenADAS
