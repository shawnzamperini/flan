/**
* @file save_results.h
* @brief Header file for save_results.cpp
*/

#ifndef SAVE_RESULTS_H
#define SAVE_RESULTS_H

#include <string>

#include "background.h"
#include "boundary.h"
#include "collisions.h"
#include "impurity_stats.h"
#include "netcdf"
#include "options.h"
#include "sputtering.h"

namespace SaveResults
{

	void save_results(const Background::Background& bkg,
		Impurity::Statistics& imp_stats, const Options::Options& opts,
		const Sputtering::Source* source = nullptr,
		const Sputtering::Deposition* deposition = nullptr,
		const Boundary::Fates* fates = nullptr,
		const std::array<Collisions::SkewCounts, 2>* skew_counts = nullptr);

	void save_netcdf(const Background::Background& bkg, 
		Impurity::Statistics& imp_stats, const Options::Options& opts,
		const Sputtering::Source* source = nullptr,
		const Sputtering::Deposition* deposition = nullptr,
		const Boundary::Fates* fates = nullptr,
		const std::array<Collisions::SkewCounts, 2>* skew_counts = nullptr);
}

#endif
