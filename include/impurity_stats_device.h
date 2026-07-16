#pragma once

namespace ImpurityStats
{
	// Device-side structure to hold data from host-side Statistics object
	struct StatisticsDevice 
	{
		int device_id;
		int Nt;
		int Nx;
		int Ny;
		int Nz;

		int *counts;
		double *weights;

		double *vX;
		double *vY;
		double *vZ;

		double *vx;
		double *vy;
		double *vz;

		double *q;
	};
}
