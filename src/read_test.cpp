/**
* @file read_test.cpp
*
* @brief Routines handling generating a test background plasma
*/
#include <cmath>

#include "background.h"
#include "flan_types.h"
#include "options.h"
#include "utilities.h"
#include "vectors.h"

namespace Test
{
	// Entry point for creating hardcoded test data for Flan
	Background::Background read_test(const Options::Options& opts)
	{	
		// The test grid is 64x64x1. 		
		constexpr int tdim {128};
		constexpr int xdim {32};
		constexpr int ydim {32};
		constexpr int zdim {9};

		// Vectors to hold the grid coordinates
		std::vector<double> test_grid_x = {};
		std::vector<double> test_grid_y = {};
		std::vector<double> test_grid_z = {};

		// Vectors to hold the cell center coordinates
		auto test_t = Utilities::linspace(0.0, 5e-6, tdim);
		std::vector<double> test_x = {};
		std::vector<double> test_y = {};
		std::vector<double> test_z = {};

		// Empty Vector3Ds to be manually filled out
		Vectors::Vector3D<BkgFPType> test_dxdX {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> test_dxdY {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> test_dxdZ {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> test_dydX {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> test_dydY {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> test_dydZ {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> test_dzdX {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> test_dzdY {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> test_dzdZ {xdim, ydim, zdim};
		Vectors::Vector3D<BkgFPType> test_J {xdim, ydim, zdim};

		// Empty Vector4Ds to be manually filled out
		//Vectors::Vector4D<BkgFPType> test_vp {tdim, xdim, ydim, zdim};
		Vectors::Vector4D<BkgFPType> test_eX {tdim, xdim, ydim, zdim};
		Vectors::Vector4D<BkgFPType> test_eY {tdim, xdim, ydim, zdim};
		Vectors::Vector4D<BkgFPType> test_eZ {tdim, xdim, ydim, zdim};
		Vectors::Vector4D<BkgFPType> test_emag {tdim, xdim, ydim, zdim};
		Vectors::Vector4D<BkgFPType> test_bX {tdim, xdim, ydim, zdim};
		Vectors::Vector4D<BkgFPType> test_bY {tdim, xdim, ydim, zdim};
		Vectors::Vector4D<BkgFPType> test_bZ {tdim, xdim, ydim, zdim};
		Vectors::Vector4D<BkgFPType> test_bmag {tdim, xdim, ydim, zdim};
		Vectors::Vector4D<BkgFPType> test_ne {tdim, xdim, ydim, zdim};
		Vectors::Vector4D<BkgFPType> test_te {tdim, xdim, ydim, zdim};
		Vectors::Vector4D<BkgFPType> test_ti {tdim, xdim, ydim, zdim};
		Vectors::Vector4D<BkgFPType> test_uX {tdim, xdim, ydim, zdim};
		Vectors::Vector4D<BkgFPType> test_uY {tdim, xdim, ydim, zdim};
		Vectors::Vector4D<BkgFPType> test_uZ {tdim, xdim, ydim, zdim};
		
		// Integers that correspond to each test
		// 0 = Simple gyration test to show a particle gyrates in a constant
		//     magnetic field
		// 1 = Test to show a gyrating particle will ExB drift in a constant
		//     electric and magnetic field
		// 2 = Test to show a gyrating particle will grad-B drift in a magentic
		//     field with a linear gradient
		// 3 = Test to show a gyrating particle will have a polarization drift
		//     in a time-varying electric field and constant magnetic field
		// 4 = Test to show a gyrating particle will have a curvature drift
		//     in a toroidal magnetic field
		// 5 = Test to show that the fluid friction force manifests  from the
		//     collision model

		// Slab geometry options (includes tests for simple gyration, ExB, 
		// grad-B, polarization drifts and friction force).
		if (opts.test_opt_int() == 0 || opts.test_opt_int() == 1 
			|| opts.test_opt_int() == 2 || opts.test_opt_int() == 3
			|| opts.test_opt_int() == 5)
		{
			// x, y, z just a simple rectangular volume, 5 x 5 x 3 cm.
			test_grid_x = Utilities::linspace(0.00, 0.05, xdim+1);
			test_grid_y = Utilities::linspace(-0.025, 0.025, ydim+1);
			test_grid_z = Utilities::linspace(-0.015, 0.015, zdim+1);

			// The cell center values are equally spaced starting at one half 
			// the cell width
			const double xwidth {test_grid_x[1] - test_grid_x[0]};
			const double ywidth {test_grid_y[1] - test_grid_y[0]};
			const double zwidth {test_grid_z[1] - test_grid_z[0]};
			test_x = Utilities::linspace(test_grid_x.front() + xwidth / 2.0, 
				test_grid_x.back() - xwidth / 2.0, xdim);	
			test_y = Utilities::linspace(test_grid_y.front() + ywidth / 2.0, 
				test_grid_y.back() - ywidth / 2.0, ydim);	
			test_z = Utilities::linspace(test_grid_z.front() + zwidth / 2.0, 
				test_grid_z.back() - zwidth / 2.0, zdim);	

			#pragma omp parallel for
			for (int i = 0; i < xdim; ++i)  // x
			{
			for (int j {}; j < ydim; ++j)  // y
			{
			for (int k {}; k < zdim; ++k)  // z
			{
				// Calculate index in 1D vector
				int idx {test_dxdX.calc_index(i,j,k)};

				// Manually assign the reciprocal basis vectors for slab, 
				// pretty straightforward since x=X, y=Y and z=Z.
				test_dxdX.get_data()[idx] = 1.0;
				test_dxdY.get_data()[idx] = 0.0;
				test_dxdZ.get_data()[idx] = 0.0;
				test_dydX.get_data()[idx] = 0.0;
				test_dydY.get_data()[idx] = 1.0;
				test_dydZ.get_data()[idx] = 0.0;
				test_dzdX.get_data()[idx] = 0.0;
				test_dzdY.get_data()[idx] = 0.0;
				test_dzdZ.get_data()[idx] = 1.0;

				// Jacobian just 1.0 since x=X, y=Y and z=Z.
				test_J.get_data()[idx] = 1.0;
			}
			}
			}
		}
	
		// Cylindrical geometry options (curvature drift and numerical
		// diffusion tests)
		else if (opts.test_opt_int() == 4)
		{
			// Normal cylindrical coordinates: x = R, y = Z, z = phi. 
			// Dimension are x = [2.00, 2.05], y = [-0.025, 0.025] z = [0, pi/2]
			test_grid_x = Utilities::linspace(2.00, 2.05, xdim+1);
			test_grid_y = Utilities::linspace(-0.025, -0.025, ydim+1);
			test_grid_z = Utilities::linspace(0.00, 3.1415 / 2.0, zdim+1);

			// The cell center values are equally spaced starting at one half 
			// the cell width
			const double xwidth {test_grid_x[1] - test_grid_x[0]};
			const double ywidth {test_grid_y[1] - test_grid_y[0]};
			const double zwidth {test_grid_z[1] - test_grid_z[0]};
			test_x = Utilities::linspace(xwidth / 2.0, 
				test_grid_x.back() - xwidth / 2.0, xdim);	
			test_y = Utilities::linspace(ywidth / 2.0, 
				test_grid_y.back() - ywidth / 2.0, ydim);	
			test_z = Utilities::linspace(zwidth / 2.0, 
				test_grid_z.back() - zwidth / 2.0, zdim);	

			#pragma omp parallel for
			for (int i = 0; i < xdim; ++i)  // x
			{
			for (int j {}; j < ydim; ++j)  // y
			{
			for (int k {}; k < zdim; ++k)  // z
			{
				// Calculate index in 1D vector
				int idx {test_dxdX.calc_index(i,j,k)};

				double R  = test_x[i];   // x = R
				//double Zc = test_y[j];   // y = Z
				double phi = test_z[k];  // z = phi
				test_dxdX.get_data()[idx] = std::cos(phi);
				test_dxdY.get_data()[idx] = std::sin(phi);
				test_dxdZ.get_data()[idx] = 0.0;
				test_dydX.get_data()[idx] = 0.0;
				test_dydY.get_data()[idx] = 0.0;
				test_dydZ.get_data()[idx] = 1.0;
				test_dzdX.get_data()[idx] = -std::sin(phi) / R;
				test_dzdY.get_data()[idx] =  std::cos(phi) / R;
				test_dzdZ.get_data()[idx] = 0.0;

				// Jacobian in cylindrical coordinates is just R
				test_J.get_data()[idx] = R;
			}
			}
			}
		}

		// Fill in time-dependent data
		#pragma omp parallel for
		for (int i = 0; i < tdim; ++i)  // t
		{
		for (int j {}; j < xdim; ++j)  // x
		{
		for (int k {}; k < ydim; ++k)  // y
		{
		for (int l {}; l < zdim; ++l)  // z
		{

			// Calculate index in 1D vector
			int idx {test_ne.calc_index(i,j,k,l)};

			// All test cases use constant values for density and temperature
			test_ne.get_data()[idx] = 1e20;
			test_te.get_data()[idx] = 1;
			test_ti.get_data()[idx] = 1;
			//test_vp.get_data()[idx] = 0.0; // Not needed

			// Zero electric field in these tests
			if (opts.test_opt_int() == 0 || opts.test_opt_int() == 2 ||
				opts.test_opt_int() == 4 || opts.test_opt_int() == 5)
			{
				test_eX.get_data()[idx] = 0.0;
				test_eY.get_data()[idx] = 0.0;
				test_eZ.get_data()[idx] = 0.0;
				test_emag.get_data()[idx] = 0.0;
			}

			// Constant electric field of E_Y = 5000 V/m for ExB drift test
			else if (opts.test_opt_int() == 1)
			{
				test_eY.get_data()[idx] = 5000.0;
				test_emag.get_data()[idx] = 5000.0;
			}

			// Time-varying electric field that starts at E_Y = 5000 V/m and
			// decreases by 1 V every microsecond for polarization drift test
			else if (opts.test_opt_int() == 3)
			{
				constexpr double slope {-1000000}; // 1000000 V/s or 1 V/us
				test_eY.get_data()[idx] = 5000.0 + slope 
					* test_t[i];
				test_emag.get_data()[idx] = test_eY.get_data()[idx];
			}

			// Constant magnetic field of B_Z = 1.0 T for gyration, ExB,
			// polarization tests.
			if (opts.test_opt_int() == 0 || opts.test_opt_int() == 1 ||
				opts.test_opt_int() == 3)
			{
				test_bX.get_data()[idx] = 0.0;
				test_bY.get_data()[idx] = 0.0;
				test_bZ.get_data()[idx] = 1.0;
				test_bmag.get_data()[idx] = 1.0;
			}

			// B_Z-gradient of 5 T/m in the Y direction for grad-B drift test
			else if (opts.test_opt_int() == 2)
			{
				constexpr double slope {5.0}; // 5 T/m or 0.05 T/cm
				test_bZ.get_data()[idx] = slope * test_y[k] + 1.0;
				test_bmag.get_data()[idx] = test_bZ.get_data()[idx];
			}

			// Cylindrical drift test gets a purely toroidal 1 T magnetic field
			else if (opts.test_opt_int() == 4)
			{
				double phi = test_z[l]; // z = phi
				double B0 = 1.0;
				test_bX.get_data()[idx] = -B0 * std::sin(phi);
				test_bY.get_data()[idx] =  B0 * std::cos(phi);
				test_bZ.get_data()[idx] =  0.0;
				test_bmag.get_data()[idx] = B0;
			}
			
			// All drift tests have zero background flow
			if (opts.test_opt_int() == 0 || opts.test_opt_int() == 1 
				|| opts.test_opt_int() == 2 || opts.test_opt_int() == 3
				|| opts.test_opt_int() == 4)
			{
				test_uX.get_data()[idx] = 0.0;
				test_uY.get_data()[idx] = 0.0;
				test_uZ.get_data()[idx] = 0.0;
			}

			// Friction force has flow of 1,000 m/s in X direction
			else if (opts.test_opt_int() == 5)
			{
				test_uX.get_data()[idx] = 1000.0;
				test_uY.get_data()[idx] = 0.0;
				test_uZ.get_data()[idx] = 0.0;
			}
		}
		}
		}
		}

		// With all our arrays assembled, encapsulate them into a Background
		// class object with move semantics then return.
		Background::Background bkg {};

		// Move each of the vectors we've created into our bkg object
		bkg.move_into_times(test_t);
		bkg.move_into_x(test_x);
		bkg.move_into_y(test_y);
		bkg.move_into_z(test_z);
		bkg.move_into_grid_x(test_grid_x);
		bkg.move_into_grid_y(test_grid_y);
		bkg.move_into_grid_z(test_grid_z);
		bkg.move_into_ne(test_ne);
		bkg.move_into_te(test_te);
		bkg.move_into_ti(test_ti);
		//bkg.move_into_vp(test_vp);
		bkg.move_into_bX(test_bX);
		bkg.move_into_bY(test_bY);
		bkg.move_into_bZ(test_bZ);
		bkg.move_into_bmag(test_bmag);
		bkg.move_into_eX(test_eX);
		bkg.move_into_eY(test_eY);
		bkg.move_into_eZ(test_eZ);
		bkg.move_into_emag(test_emag);
		bkg.move_into_uX(test_uX);
		bkg.move_into_uY(test_uY);
		bkg.move_into_uZ(test_uZ);

		// Likewise for the reciprocal basis vectors. Tangent basis vectors
		// not needed, so not saved.
		bkg.move_into_dxdX(test_dxdX);
		bkg.move_into_dxdY(test_dxdY);
		bkg.move_into_dxdZ(test_dxdZ);
		bkg.move_into_dydX(test_dydX);
		bkg.move_into_dydY(test_dydY);
		bkg.move_into_dydZ(test_dydZ);
		bkg.move_into_dzdX(test_dzdX);
		bkg.move_into_dzdY(test_dzdY);
		bkg.move_into_dzdZ(test_dzdZ);
		bkg.move_into_J(test_J);

		// Tangent basis vectors
		//bkg.move_into_dXdx(test_dXdx);
		//bkg.move_into_dYdx(test_dYdx);
		//bkg.move_into_dZdx(test_dZdx);
		//bkg.move_into_dXdy(test_dXdy);
		//bkg.move_into_dYdy(test_dYdy);
		//bkg.move_into_dZdy(test_dZdy);
		//bkg.move_into_dXdz(test_dXdz);
		//bkg.move_into_dYdz(test_dYdz);
		//bkg.move_into_dZdz(test_dZdz);

		// Okay to return by value since C++11 uses move semantics here.
		return bkg;
	}
}
