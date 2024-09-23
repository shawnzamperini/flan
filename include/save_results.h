/**
* @file save_results.h
* @brief save_results.cpp header file
*/

#ifndef SAVE_RESULTS_H
#define SAVE_RESULTS_H

#include <string>

#include "background.h"
#include "impurity_stats.h"
#include "netcdf"

namespace SaveResults
{

	// Entry point to saving results.
	void save_results(const std::string_view case_name, 
		const Background::Background& bkg, Impurity::Statistics& imp_stats);

	
	// Save a 1D vector into nc_file. The correct NcDim (dim) must be passed 
	// in - it on the programmer to do this correctly.
	template <typename T>
	void save_vector_1d(const netCDF::NcFile& nc_file, std::vector<T> vec, 
		const std::string& var_name, const netCDF::NcDim& dim, 
		const std::string& description, const std::string& units);

	
	// Save a 4D vector into nc_file. The correct NcDim (dim) must be passed 
	// in - it on the programmer to do this correctly. Currently all 
	// Vector4D's hold double, so no template like in save_vector_1d.
	template <typename T>
	void save_vector_4d(const netCDF::NcFile& nc_file, 
		const Vectors::Vector4D<T>& vec, 
		const std::string& var_name, const netCDF::NcDim& dim1,
		const netCDF::NcDim& dim2, const netCDF::NcDim& dim3, 
		const netCDF::NcDim& dim4, const std::string& description, 
		const std::string& units);

	// Save all results in a netCDF4 file.
	void save_netcdf(const std::string_view case_name, 
		const Background::Background& bkg, Impurity::Statistics& imp_stats);

}

#endif
