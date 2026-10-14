/**
* @file openadas.h
* @brief Header file for openadas.cpp
*/
#pragma once

#include <string>
#include <vector>

#include "background.h"
#include "openadas_device.h"
#include "pcg32.h"
#include "slots.h"
#include "vectors.h"


namespace OpenADAS
{

	/**
	* Class to read and store the ionization and recombination rate
	* coefficients from OpenADAS.
	*/
	class OpenADAS
	{
	private:
		int m_atomic_number {};
		int m_ndens {};
		int m_ntemp {};
		int m_charge_low {};
		int m_charge_high {};
		std::vector<double> m_te {};
		std::vector<double> m_ne {};
		Vectors::Vector3D<double> m_rates {};

	public:
		
		OpenADAS();

		OpenADAS(const std::string_view openadas_root, const int openadas_year, 
			const int imp_atom_num, const std::string_view rate_type);
		
		std::string get_openadas_path(const std::string_view openadas_root, 
			const int openadas_year, const int imp_atom_num, 
			const std::string_view rate_type);

		void read_rate_coefficients(std::ifstream& openadas_stream);

		std::tuple<double, double, double, int, int> 
			get_bilinear_interp_vals(const std::vector<double>& vec, 
			double value) const;

		double get_rate_coeff(const int charge, double ne, double te) const;


		// Create a device-side struct with the OpenADAS data. Defined in
		// src/cuda/openadas.cu
		OpenADASDevice to_device(int device_id);
	};


	void ioniz_recomb(Slots::Slots& slots, 
		const Background::Background& bkg, const OpenADAS& oa_ioniz, 
		const OpenADAS& oa_recomb, const double dt, 
		int& ioniz_warnings, int& recomb_warnings, std::vector<pcg32>& rngs);

	// Free memory of device-side OpenADAS struct
	void free_oa(OpenADASDevice& oa_d, int device_id);
}
