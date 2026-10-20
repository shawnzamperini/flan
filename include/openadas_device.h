/**
* @file  openadas_device.h
* @brief Contains structs to hold OpenADAS rate coefficients and the
*        warning counters for device-side access.
*/
#pragma once

namespace OpenADAS
{
	/**
	* Device-side copy of one OpenADAS rate file (ionization or
	* recombination), made by OpenADAS::to_device and freed with free_oa.
	*
	* The values are already converted from the file's log10 units. rates
	* holds one ntemp x ndens table per charge state, one after the other,
	* each indexed by te_index * ndens + ne_index (see get_charge_slice in
	* openadas.cu).
	*/
	struct OpenADASDevice
	{
		int device_id;		///< ID of the GPU the arrays are on

		int atomic_number;	///< Atomic number of the element
		int ndens;			///< Number of densities in the table
		int ntemp;			///< Number of temperatures in the table
		int charge_low;		///< Lowest charge in the file (IZ1MIN)
		int charge_high;	///< Highest charge in the file (IZ1MAX)

		// Device pointers
		double* te;			///< Electron temperatures (eV), size ntemp
		double* ne;			///< Electron densities (m^-3), size ndens
		double* rates;		///< Rate coefficients (m^3/s), size ndens *
							///< ntemp * (charge_high - charge_low + 1)

	};  // OpenADASDevice

	/**
	* Device-side counters of how many times the ionization or
	* recombination probability exceeded 1, and of how many particles were
	* killed for exceeding imp_max_time.
	*
	* Allocated once per GPU with alloc_warnings, incremented by
	* ioniz_recomb_kernel (ioniz, recomb) and ImpurityTransport::step_kernel
	* (max_time) every time step, read back with add_warnings_to_host at
	* the end of the run and freed with free_warnings.
	*/
	struct IonizRecombWarningsDevice
	{
		int device_id {};	///< ID of the GPU the counters are on
		int* ioniz {};		///< Single int, number of times the
							///< ionization probability exceeded 1
		int* recomb {};		///< Single int, number of times the
							///< recombination probability exceeded 1
		int* max_time {};	///< Single int, number of particles killed
							///< for exceeding imp_max_time
	};  // IonizRecombWarningsDevice

}  // namespace OpenADAS
