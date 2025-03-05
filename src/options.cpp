/**
* @file options.cpp
*/
#include "options.h"

namespace Options
{
	// Constructors
	Options::Options()
	{}
	
	// Setter definitions
	void Options::set_imp_num(int imp_num) {m_imp_num = imp_num;}

	// Accessor definitions
	int Options::imp_num() {return m_imp_num;}
}
