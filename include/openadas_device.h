#pragma once

namespace OpenADAS
{
	struct OpenADASDevice
	{
		int device_id;

		int atomic_number;
		int ndens;
		int ntemp;
		int charge_low;
		int charge_high;

		// Device pointers
		double* te;      // size ntemp
		double* ne;      // size ndens
		double* rates;   // size ndens * ntemp * (charge_high - charge_low + 1)

	};  // OpenADASDevice

	/**
	* Device-side counters of how many times the ionization or
	* recombination probability exceeded 1.
	*
	* Allocated once per GPU with alloc_warnings, incremented by
	* ioniz_recomb_kernel every time step, read back with
	* add_warnings_to_host at the end of the run and freed with
	* free_warnings.
	*/
	struct IonizRecombWarningsDevice
	{
		int device_id {};
		int* ioniz {};    // single int on the device
		int* recomb {};   // single int on the device
	};  // IonizRecombWarningsDevice

}  // namespace OpenADAS
