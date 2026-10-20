#pragma once

/**
* @file  read_bkg.h
*
* @brief Header file for read_bkg.cpp
*/
#include "background.h"
#include "options.h"

namespace Background
{
	Background read_bkg(const Options::Options& opts);
}
