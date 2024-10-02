/**
* @file openadas.cpp
*/
#include <vector>
#include <string>
#include <iostream>
#include <cassert>

#include "openadas.h"
#include "vectors.h"


namespace OpenADAS
{

	OpenADAS::OpenADAS(std::string_view openadas_root, int openadas_year, 
		int imp_atom_num, std::string_view rate_type)
	{
		// Create full path to recombination or ionization rate files. 
		// These are the ACD and SCD files, respectively.
		std::string openadas_path {get_openadas_path(openadas_root, 
			openadas_year, imp_atom_num, rate_type)};

		// Open file

		// Read in electron temperatures and densities that data is 
		// availiable for into vectors. 

		// Read in rate coefficients into Vector3D.
	}

	std::string OpenADAS::get_openadas_path(std::string_view openadas_root, 
		int openadas_year, int imp_atom_num, std::string_view rate_type)
	{
		// Check for valid entry for rate type
		if (rate_type != "scd" && rate_type != "acd")
		{
			std::cerr << "Error! rate_type = " << rate_type << ". Only " 
				<< "\"scd\" and \"acd\" are allowed.\n";
			// throw exception...
		}

		// Select element name from atomic number.
		std::string imp_name {};
		if (imp_atom_num == 12) imp_name = "c";
		else if (imp_atom_num == 74) imp_name = "w";
		else
		{
			std::cerr << "Error! imp_atom_num = " << imp_atom_num << " is not"
				<< " a valid option.\n";
			// throw exception...
		}

		// Assemble the full path. Example: /home/zamp/openadas/scd50_w.dat
		std::string openadas_path {openadas_root};
		openadas_path.append("/");
		openadas_path.append(rate_type);
		openadas_path.append(std::to_string(openadas_year));
		openadas_path.append("_");
		openadas_path.append(imp_name); 
		openadas_path.append(".dat");

		return openadas_path;
	}
}

