#include <iostream>
#include <tuple>

#include "flan.h"
#include "flan_types.h"


// Function definitions go here (don't forget to make a pointer below)
std::tuple<double, double, double> mapc2p(double xc, double yc, double zc)
{
	// Default template where the computational coordinates are already
	// in physical (Cartesian) coordinates. Modify as needed based on the
	// Gkeyll simulation.
	double xp {xc};
	double yp {yc};
	double zp {zc};

	std::cout << "Using passed in mapc2p\n";

	// Return as three-tuple of doubles.
	return std::make_tuple(xp, yp, zp);
}

Inputs create_inputs()
{
	Inputs inpts;
	
	// Int, double and string input options go here
	inpts["imp_num"] = 10;

	// Pointers to functions go here
	Mapc2p_ptr mapc2p_ptr {&mapc2p};

	// Function pointer input options go here
	inpts["mapc2p"] = mapc2p_ptr;

	return inpts;
}


int main()
{
	// Run Flan
	Inputs inpts {create_inputs()};
	flan(inpts);
}
