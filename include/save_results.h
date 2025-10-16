/**
* @file save_results.h
* @brief Header file for save_results.cpp
*/

#ifndef SAVE_RESULTS_H
#define SAVE_RESULTS_H

#include <string>

#include "background.h"
#include "impurity_stats.h"
#include "netcdf"
#include "options.h"

namespace SaveResults
{

	// Entry point to saving results.
	void save_results(const Background::Background& bkg,
		Impurity::Statistics& imp_stats, const Options::Options& opts);

	// Save all results in a netCDF4 file.
	void save_netcdf(const Background::Background& bkg, 
		Impurity::Statistics& imp_stats, const Options::Options& opts);
}

#endif
