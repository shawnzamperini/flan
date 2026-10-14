#pragma once

#include "background.h"
#include "options.h"

namespace Background
{
	Background read_bkg(const Options::Options& opts);
}
