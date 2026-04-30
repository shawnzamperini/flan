/**
* @file read_test.cpp
*
* @brief Routines handling generating a test background plasma
*/

#include "background.h"
#include "flan_types.h"
#include "options.h"
#include "vectors.h"

namespace Test
{
	// Entry point for creating hardcoded test data for Flan
	Background::Background read_test(const Options::Options& opts)
	{	
		// The test grid is 64x64x1. The dimensions of the grid depend on if
		// it is a slab or cylindrical geometry.
		constexpr int tdim {100};
		constexpr int xdim {64};
		constexpr int ydim {64};
		constexpr int zdim {1};

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

		// Empty Vector4Ds to be manually filled out
		Vectors::Vector4D<BkgFPType> test_ne {tdim, xdim, ydim, zdim};
		Vectors::Vector4D<BkgFPType> test_te {tdim, xdim, ydim, zdim};
		Vectors::Vector4D<BkgFPType> test_ti {tdim, xdim, ydim, zdim};
		

		// Set spatial, non-time-dependent data
		#pragma omp parallel for
		for (int i = 0; i < xdim; ++i)  // x
		{
		for (int j {}; j < ydim; ++j)  // y
		{
		for (int k {}; k < zdim; ++k)  // z
		{
			// Calculate index in 1D vector
			int idx {test_dxdX.calc_index(i,j,k)};

			// Need to manually assign the reciprocal basis vectors to resemble
			// a specific geometry

			// Slab
			test_dxdX.get_data()[idx] = 1.0;
			test_dxdY.get_data()[idx] = 0.0;
			test_dxdZ.get_data()[idx] = 0.0;
			test_dydX.get_data()[idx] = 0.0;
			test_dydY.get_data()[idx] = 1.0;
			test_dydZ.get_data()[idx] = 0.0;
			test_dzdX.get_data()[idx] = 0.0;
			test_dzdY.get_data()[idx] = 0.0;
			test_dzdZ.get_data()[idx] = 1.0;

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

			// Overwrite with constant plasma values
			test_ne.get_data()[idx] = 1e20;
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

		// With all our arrays assembled, encapsulate them into a Background
		// class object with move semantics then return.
		Background::Background bkg {};

		// Move each of the vectors we've created into our bkg object
		//bkg.move_into_times(gkyl_times);
		//bkg.move_into_x(gkyl_x);
		//bkg.move_into_y(gkyl_y);
		//bkg.move_into_z(gkyl_z);
		//bkg.move_into_grid_x(gkyl_grid_x);
		//bkg.move_into_grid_y(gkyl_grid_y);
		//bkg.move_into_grid_z(gkyl_grid_z);
		bkg.move_into_ne(test_ne);
		//bkg.move_into_te(gkyl_te);
		//bkg.move_into_ti(gkyl_ti);
		//bkg.move_into_vp(gkyl_vp);
		//bkg.move_into_bX(gkyl_bX);
		//bkg.move_into_bY(gkyl_bY);
		//bkg.move_into_bZ(gkyl_bZ);
		//bkg.move_into_bmag(gkyl_bmag);
		//bkg.move_into_eX(gkyl_eX);
		//bkg.move_into_eY(gkyl_eY);
		//bkg.move_into_eZ(gkyl_eZ);
		//bkg.move_into_emag(gkyl_emag);
		//bkg.move_into_uX(gkyl_uX);
		//bkg.move_into_uY(gkyl_uY);
		//bkg.move_into_uZ(gkyl_uZ);
		//bkg.move_into_X(gkyl_X);
		//bkg.move_into_Y(gkyl_Y);
		//bkg.move_into_Z(gkyl_Z);
		//bkg.move_into_grid_X(gkyl_grid_X);
		//bkg.move_into_grid_Y(gkyl_grid_Y);
		//bkg.move_into_grid_Z(gkyl_grid_Z);
		
		// Only needed for debuggin purposes, delete if not debugging
		//bkg.move_into_gradbX(gkyl_dbdX);
		//bkg.move_into_gradbY(gkyl_dbdY);
		//bkg.move_into_gradbZ(gkyl_dbdZ);

		// Special case. gkyl_J is a 4D vector so we didn't have to write a
		// whole new function just for a 3D vector - less maintainable code. 
		// So we just need to index any random timeslice as a Vector3D, pass
		// it in. 
		//Vectors::Vector3D gkyl_J_3D {gkyl_J.slice_dim1(0)};
		//bkg.move_into_J(gkyl_J_3D);

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

		// Likewise for covariant components of magnetic field
		//Vectors::Vector3D gkyl_b_x_3D {gkyl_b_x.slice_dim1(0)};
		//bkg.move_into_b_x(gkyl_b_x_3D);
		//Vectors::Vector3D gkyl_b_y_3D {gkyl_b_y.slice_dim1(0)};
		//bkg.move_into_b_y(gkyl_b_y_3D);
		//Vectors::Vector3D gkyl_b_z_3D {gkyl_b_z.slice_dim1(0)};
		//bkg.move_into_b_z(gkyl_b_z_3D);

		// Okay to return by value since C++11 uses move semantics here.
		return bkg;
	}
}
