#pragma once

/**
* @file  read_test.h
*
* @brief Header file for read_test.cpp
*/
#include "background.h"
#include "options.h"

namespace Test
{
	Background::Background read_test(const Options::Options& opts);
}
