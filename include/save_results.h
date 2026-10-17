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

	void save_results(const Background::Background& bkg,
		Impurity::Statistics& imp_stats, const Options::Options& opts);

	void save_netcdf(const Background::Background& bkg, 
		Impurity::Statistics& imp_stats, const Options::Options& opts);
}

#endif
