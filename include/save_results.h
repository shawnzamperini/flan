#ifndef SAVE_RESULTS_H
#define SAVE_RESULTS_H

#include <string>

#include "background.h"

namespace SaveResults
{

	void save_results(const std::string_view case_name, const Background::Background& bkg);

	void save_netcdf(const std::string_view case_name, const Background::Background& bkg);

}

#endif
