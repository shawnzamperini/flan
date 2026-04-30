/**
* @file read_test.cpp
*
* @brief Routines handling generating a test background plasma
*/

#include "background.h"
#include "options.h"

namespace Test
{
	// Entry point for creating hardcoded test data for Flan
	Background::Background read_test(const Options::Options& opts)
	{	

/*
	// ---- BEGIN HARDCODED TEST SUITE ----
	// Flan is missing a test suite. For the paper, the reviewer wants
	// some kind of test, so I am going to overwrite the background with 
	// constant values to test simple particle gyromotion.

	#pragma omp parallel for
	for (int i = 0; i < std::ssize(bkg.get_x()); ++i)  // x
	{
	for (int j {}; j < std::ssize(bkg.get_y()); ++j)  // y
	{
	for (int k {}; k < std::ssize(bkg.get_z()); ++k)  // z
	{
		// Calculate index in 1D vector
		//int idx {bkg.get_dxdX().calc_index(i,j,k)};

		// Need to manually assign the reciprocal basis vectors to resemble
		// a specific geometry

		// Slab
		//bkg.get_dxdX().get_data()[idx] = 1.0;
		//bkg.get_dxdY().get_data()[idx] = 0.0;
		//bkg.get_dxdZ().get_data()[idx] = 0.0;
		//bkg.get_dydX().get_data()[idx] = 0.0;
		//bkg.get_dydY().get_data()[idx] = 1.0;
		//bkg.get_dydZ().get_data()[idx] = 0.0;
		//bkg.get_dzdX().get_data()[idx] = 0.0;
		//bkg.get_dzdY().get_data()[idx] = 0.0;
		//bkg.get_dzdZ().get_data()[idx] = 1.0;

		// Sheared slab
		//bkg.get_dxdX().get_data()[idx] = 1.0;
		//bkg.get_dxdY().get_data()[idx] = 0.0;
		//bkg.get_dxdZ().get_data()[idx] = 0.0;
		//bkg.get_dydX().get_data()[idx] = 0.0;
		//bkg.get_dydY().get_data()[idx] = 1.0;
		//bkg.get_dydZ().get_data()[idx] = 0.0;
		//bkg.get_dzdX().get_data()[idx] = -3.0;   // shear
		//bkg.get_dzdY().get_data()[idx] = 0.0;
		//bkg.get_dzdZ().get_data()[idx] = 1.0;

		// Cylindrical: Reassign z to be toroidal angle from 0-pi/2 (phi).
		//bkg.get_z()[k] = static_cast<double>(k) / std::ssize(bkg.get_z()) 
			//* 3.1415 / 2.0;
		//double R  = bkg.get_x()[i];   // x = R
		//double Zc = bkg.get_y()[j];   // y = Z
		//double phi = bkg.get_z()[k];  // z = phi
		//bkg.get_dxdX().get_data()[idx] = std::cos(phi);
		//bkg.get_dxdY().get_data()[idx] = std::sin(phi);
		//bkg.get_dxdZ().get_data()[idx] = 0.0;
		//bkg.get_dydX().get_data()[idx] = 0.0;
		//bkg.get_dydY().get_data()[idx] = 0.0;
		//bkg.get_dydZ().get_data()[idx] = 1.0;
		//bkg.get_dzdX().get_data()[idx] = -std::sin(phi) / R;
		//bkg.get_dzdY().get_data()[idx] =  std::cos(phi) / R;
		//bkg.get_dzdZ().get_data()[idx] = 0.0;
	}
	}
	}

	#pragma omp parallel for
	for (int i = 0; i < std::ssize(bkg.get_times()); ++i)  // t
	{
	for (int j {}; j < std::ssize(bkg.get_x()); ++j)  // x
	{
	for (int k {}; k < std::ssize(bkg.get_y()); ++k)  // y
	{
	for (int l {}; l < std::ssize(bkg.get_z()); ++l)  // z
	{

		// Calculate index in 1D vector
		//int idx {bkg.get_ne().calc_index(i,j,k,l)};

		// Overwrite with constant plasma values
		//bkg.get_ne().get_data()[idx] = 1e20;
		//bkg.get_te().get_data()[idx] = 1;
		//bkg.get_ti().get_data()[idx] = 1;
		//bkg.get_vp().get_data()[idx] = 0.0; // Not needed

		// No electric field
		//bkg.get_eX().get_data()[idx] = 0.0;
		//bkg.get_eY().get_data()[idx] = 0.0;
		//bkg.get_eZ().get_data()[idx] = 0.0;
		//bkg.get_emag().get_data()[idx] = 0.0;

		// Slab: 500 V/m E field in Y direction
		//bkg.get_eY().get_data()[idx] = 500.0;
		//bkg.get_emag().get_data()[idx] = 500.0;

		// Slab: Time varying E field for polarization drift
		//constexpr double slope {-1000000}; // 1000000 V/s or 1 V/us
		//bkg.get_eY().get_data()[idx] = 500.0 + slope * bkg.get_times()[i];
		//bkg.get_emag().get_data()[idx] = 500.0 + slope * bkg.get_times()[i];

		// Slab: Magnetic field. Choose correct bmag.
		//bkg.get_bX().get_data()[idx] = 0.0;
		//bkg.get_bX().get_data()[idx] = 1.0;
		//bkg.get_bY().get_data()[idx] = 0.0;
		//bkg.get_bZ().get_data()[idx] = 0.0;
		//bkg.get_bmag().get_data()[idx] = 1.0;
		//bkg.get_bmag().get_data()[idx] = 0.0;
		
		// Slab: BZ-gradient in y direction for grad-B drift
		//constexpr double slope {5}; // 5 T/m or 0.05 T/cm
		//bkg.get_bZ().get_data()[idx] = slope * bkg.get_y()[k] + 1.0;
		//bkg.get_bmag().get_data()[idx] = slope * bkg.get_y()[k] + 1.0;

		// Slab: No flow or set to a value for friction force testing
		//bkg.get_uX().get_data()[idx] = 0.0;
		//bkg.get_uX().get_data()[idx] = 1000.0;
		//bkg.get_uY().get_data()[idx] = 0.0;
		//bkg.get_uZ().get_data()[idx] = 0.0;

		// Cylindrical: Purely toroidal B field
		//double phi = bkg.get_z()[l]; // z = phi
		//double B0 = 1.0;
		//bkg.get_bX().get_data()[idx] = -B0 * std::sin(phi);
		//bkg.get_bY().get_data()[idx] =  B0 * std::cos(phi);
		//bkg.get_bZ().get_data()[idx] =  0.0;
		//bkg.get_bmag().get_data()[idx] = B0;
	}
	}
	}
	}
	// ---- END HARDCODED TEST SUITE ----
*/
	}
}
