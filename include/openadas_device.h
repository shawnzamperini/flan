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

}  // namespace OpenADAS
