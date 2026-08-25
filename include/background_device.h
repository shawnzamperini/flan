#pragma once 

namespace Background
{
	struct BackgroundDevice
	{
		int device_id;
		int tdim;
		int xdim;
		int ydim;
		int zdim;

		// These are stored in constant memory, see cuda/device_constants.cuh.
		// They are in the Background namespace and can be accessed as d_t,
		// d_x, d_y, etc. 
		//double* times;
		//double* x;
		//double* y;
		//double* z;
		//double* grid_x;
		//double* grid_y;
		//double* grid_z;

		double* ne;
		double* te;
		double* ti;
		//double* vp;

		double* bX;
		double* bY;
		double* bZ;
		double* bmag;

		//double* gradbX;
		//double* gradbY;
		//double* gradbZ;

		double* eX;
		double* eY;
		double* eZ;
		double* emag;

		double* uX;
		double* uY;
		double* uZ;

		double* X;
		double* Y;
		double* Z;

		double* grid_X;
		double* grid_Y;
		double* grid_Z;

		double* J;

		double* gij_00;
		double* gij_01;
		double* gij_02;
		double* gij_11;
		double* gij_12;
		double* gij_22;

		double* dxdX;
		double* dxdY;
		double* dxdZ;

		double* dydX;
		double* dydY;
		double* dydZ;

		double* dzdX;
		double* dzdY;
		double* dzdZ;

		double* dXdx;
		double* dYdx;
		double* dZdx;

		double* dXdy;
		double* dYdy;
		double* dZdy;

		double* dXdz;
		double* dYdz;
		double* dZdz;

	};  // BackgroundDevice

}  // namespace Background
