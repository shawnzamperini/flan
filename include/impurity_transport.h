/**
* @file impurity_transport.h
*
* @brief impurity_transport.cpp header file
*/

#ifndef IMPURITY_TRANSPORT_H
#define IMPURITY_TRANSPORT_H

#include "background.h"
#include "impurity_stats.h"

namespace Impurity
{
	Statistics follow_impurities(Background::Background& bkg);
}

#endif
