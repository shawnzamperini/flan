/**
* @file openadas.h
* @brief Header file for openadas.cpp
*/
#ifndef OPENADAS_H
#define OPENADAS_H

#include <vector>
#include <string>

#include "vectors.h"
#include "impurity.h"
#include "background.h"


namespace OpenADAS
{

	/**
	* @brief Class to read and store ionization/recombination rates from 
	* OpenADAS
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
		
		/**
		* @brief Constructor
		*
		* @param openadas_root Full path to the directory containing OpenADAS 
		* data
		* @param openadas_year The year of the OpenADAS data to load
		* @param imp_atom_num Atomic number of the ion to load data for
		* @param rate_type One of "acd" (recombination) or "scd" (ionization)
		*/
		OpenADAS(const std::string_view openadas_root, const int openadas_year, 
			const int imp_atom_num, const std::string_view rate_type);
		
		/**
		* @brief Construct full path to an OpenADAS data file
		*
		* @param openadas_root Full path to the directory containing OpenADAS
		* data
		* @param openadas_year The year of the OpenADAS data to load
		* @param imp_atom_num Atomic number of the ion to load data for
		* @param rate_type One of "acd" (recombination) or "scd" (ionization)
		*
		* @return Returns string containing the full path to the OpenADAS file
		* that data is being read in from.
		*/
		std::string get_openadas_path(const std::string_view openadas_root, 
			const int openadas_year, const int imp_atom_num, 
			const std::string_view rate_type);

		/**
		* @brief Read in data from an OpenADAS file, storing it internally
		* to the class
		*
		* @param openadas_stream A filestream object of the OpenADAS file
		*/
		void read_rate_coefficients(std::ifstream& openadas_stream);

		/**
		*
		*/
		std::tuple<double, double, double, int, int> 
			get_bilinear_interp_vals(const std::vector<double>& vec, 
			double value) const;

		/**
		*
		*/
		double get_rate_coeff(const int charge, double ne, double te) const;
	};

	/**
	* @brief Calculate ionization and recombination probabilities
	* @return Returns a pair of doubles representing the probabilities of each
	*	process occurring (ionization, recombination).
	*/
	std::pair<double, double> calc_ioniz_recomb_probs(
		Impurity::Impurity& imp, const Background::Background& bkg,
		const OpenADAS& oa_ioniz, const OpenADAS& oa_recomb, 
		const double imp_time_step, const int tidx, const int xidx, 
		const int yidx, const int zidx);

	/**
	* @brief Handle impurity ion ionization and recombination
	*
	* The impurity charge is modified within this function, if necessary. If
	* the probabilities for ionization or recombination are greater than 1.0,
	* ioniz_warnings and recomb_warnings are incremented by 1. 
	*
	* @param ioniz_warnings Integer that tracks number of ionization 
	* probability > 1.0 events
	* @param recomb_warnings Integer that tracks number of recombination 
	* probability > 1.0 events
	*/
	void ioniz_recomb(Impurity::Impurity& imp, 
		const Background::Background& bkg, const OpenADAS& oa_ioniz, 
		const OpenADAS& oa_recomb, const double imp_time_step, 
		const int tidx, const int xidx, const int yidx, const int zidx, 
		int& ioniz_warnings, int& recomb_warnings);
}

#endif
