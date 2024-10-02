/**
* @file openadas.h
*/
#ifndef OPENADAS_H
#define OPENADAS_H

#include <vector>
#include <string>

#include "vectors.h"


namespace OpenADAS
{

	class OpenADAS
	{
	private:
		std::vector<double> m_te {};
		std::vector<double> m_ne {};
		Vectors::Vector3D<double> m_rates {};

	public:
		
		/**
		* @brief Constructor
		*
		* @param openadas_root Full path to the directory containing OpenADAS 
		* data
		* @param imp_atom_num Atomic number of the ion to load data for
		* @param rate_type One of "acd" (recombination) or "scd" (ionization)
		*/
		OpenADAS(std::string_view openadas_root, int openadas_year, 
			int imp_atom_num, std::string_view rate_type);
		
		/**
		*
		*/
		std::string get_openadas_path(std::string_view openadas_root, 
			int openadas_year, int imp_atom_num, std::string_view rate_type);
	};
}

#endif
