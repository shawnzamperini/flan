#pragma once

/**
* @file  read_bkg.h
*
* @brief Header file for read_bkg.cpp
*/
#include <vector>

#include "background.h"
#include "flan_types.h"
#include "options.h"
#include "zones.h"

namespace Background
{
	Background read_bkg(const Options::Options& opts);

	Zones::Zones single_zone(std::vector<BkgFPType>& x,
		std::vector<BkgFPType>& y, std::vector<BkgFPType>& z,
		std::vector<BkgFPType>& grid_x, std::vector<BkgFPType>& grid_y,
		std::vector<BkgFPType>& grid_z, const Options::Options& opts);
}
