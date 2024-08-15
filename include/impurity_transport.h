#ifndef IMPURITY_TRANSPORT_H
#define IMPURITY_TRANSPORT_H

#include "read_gkyl.h"  // For Background

namespace Impurity
{
	void follow_impurities(Gkyl::Background& bkg);
}

#endif
