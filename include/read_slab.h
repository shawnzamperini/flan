#pragma once

/**
* @file  read_slab.h
*
* @brief Header file for read_slab.cpp
*/
#include "background.h"
#include "options.h"

namespace Slab
{
	Background::Background read_slab(const Options::Options& opts);
}
