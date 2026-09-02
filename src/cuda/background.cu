#include <cuda_runtime.h> 

#include "background.h"
#include "background_device.h"

#include "device_constants.cuh"


namespace Background
{

	// Copy background data to device, returning a BackgroundDevice object with
	// pointers to memory location on device.
	BackgroundDevice Background::to_device(int device_id)
	{
		BackgroundDevice bkg_d {};

	#ifdef USE_CUDA

		bkg_d.device_id = device_id;

		// In include/cuda/device_constants.cuh we define arrays like t, x, 
		// etc. in constant memory since they're small. Their size is set by
		// MAX_1D in that file. As time goes on this will need to be increased,
		// but for now print a error that if any dimension is above that you
		// WILL run into difficult to find errors due to accessing out of
		// bounds memory (which Nvidia GPUs don't care about).
		if (m_dim1 > MAX_1D || m_dim2 > MAX_1D
			|| m_dim3 > MAX_1D || m_dim4 > MAX_1D)
		{
			std::cerr << "Error! MAX_1D is too small. It MUST be increased.\n"
				<< "  dim1 = " << m_dim1 << "\n"
				<< "  dim2 = " << m_dim2 << "\n"
				<< "  dim3 = " << m_dim3 << "\n"
				<< "  dim4 = " << m_dim4 << "\n";
		}

		bkg_d.tdim = m_dim1;
		bkg_d.xdim = m_dim2;
		bkg_d.ydim = m_dim3;
		bkg_d.zdim = m_dim4;

		int N_3D {m_dim2 * m_dim3 * m_dim4};
		int N_4D {m_dim1 * m_dim2 * m_dim3 * m_dim4};

		cudaSetDevice(device_id);

		// -----------------------------------
		// Constant memory copies
		// -----------------------------------

		// These arrays are great candidates for putting in GPU constant
		// memory since they are small. Arrays in constant memory can be read
		// significantly faster, which is convienent since these are read
		// very often. We declare them as constant in in 
		// cuda/device_constants.cuh, and use cudaMemcpyToSymbol to copy them 
		// into constant memory.
		cudaMemcpyToSymbol(d_t, m_times.data(), m_dim1 * sizeof(double));
		cudaMemcpyToSymbol(d_x, m_x.data(), m_dim2 * sizeof(double));
		cudaMemcpyToSymbol(d_y, m_y.data(), m_dim3 * sizeof(double));
		cudaMemcpyToSymbol(d_z, m_z.data(), m_dim4 * sizeof(double));
		cudaMemcpyToSymbol(d_grid_x, m_grid_x.data(), (m_dim2+1) 
			* sizeof(double));
		cudaMemcpyToSymbol(d_grid_y, m_grid_y.data(), (m_dim3+1) 
			* sizeof(double));
		cudaMemcpyToSymbol(d_grid_z, m_grid_z.data(), (m_dim4+1) 
			* sizeof(double));

		// Copying scalars to constant memory is a little tricker since they
		// are not treated the same way as arrays. So we use this overload
		// of cudaMemcpyToSymbol which handles it.
		const double t_min {m_times.front()};
		const double t_max {m_times.back()};
		const double x_min {m_x.front()};
		const double x_max {m_x.back()};
		const double y_min {m_y.front()};
		const double y_max {m_y.back()};
		const double z_min {m_z.front()};
		const double z_max {m_z.back()};
		cudaMemcpyToSymbol(d_t_min, &t_min, sizeof(double));
		cudaMemcpyToSymbol(d_t_max, &t_max, sizeof(double));
		cudaMemcpyToSymbol(d_x_min, &x_min, sizeof(double));
		cudaMemcpyToSymbol(d_x_max, &x_max, sizeof(double));
		cudaMemcpyToSymbol(d_y_min, &y_min, sizeof(double));
		cudaMemcpyToSymbol(d_y_max, &y_max, sizeof(double));
		cudaMemcpyToSymbol(d_z_min, &z_min, sizeof(double));
		cudaMemcpyToSymbol(d_z_max, &z_max, sizeof(double));

		// -----------------------------------
		// Global memory allocation and copies
		// -----------------------------------

		cudaMalloc(&bkg_d.dxdX, N_3D * sizeof(double));
		cudaMalloc(&bkg_d.dxdY, N_3D * sizeof(double));
		cudaMalloc(&bkg_d.dxdZ, N_3D * sizeof(double));
		cudaMemcpy(bkg_d.dxdX, m_dxdX.get_data().data(), N_3D 
			* sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(bkg_d.dxdY, m_dxdY.get_data().data(), N_3D 
			* sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(bkg_d.dxdZ, m_dxdZ.get_data().data(), N_3D 
			* sizeof(double), cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.dydX, N_3D * sizeof(double));
		cudaMalloc(&bkg_d.dydY, N_3D * sizeof(double));
		cudaMalloc(&bkg_d.dydZ, N_3D * sizeof(double));
		cudaMemcpy(bkg_d.dydX, m_dydX.get_data().data(), N_3D 
			* sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(bkg_d.dydY, m_dydY.get_data().data(), N_3D 
			* sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(bkg_d.dydZ, m_dydZ.get_data().data(), N_3D 
			* sizeof(double), cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.dzdX, N_3D * sizeof(double));
		cudaMalloc(&bkg_d.dzdY, N_3D * sizeof(double));
		cudaMalloc(&bkg_d.dzdZ, N_3D * sizeof(double));
		cudaMemcpy(bkg_d.dzdX, m_dzdX.get_data().data(), N_3D 
			* sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(bkg_d.dzdY, m_dzdY.get_data().data(), N_3D 
			* sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(bkg_d.dzdZ, m_dzdZ.get_data().data(), N_3D 
			* sizeof(double), cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.ne, N_4D * sizeof(double));
		cudaMemcpy(bkg_d.ne, m_ne.get_data().data(), N_4D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.te, N_4D * sizeof(double));
		cudaMemcpy(bkg_d.te, m_te.get_data().data(), N_4D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.ti, N_4D * sizeof(double));
		cudaMemcpy(bkg_d.ti, m_ti.get_data().data(), N_4D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.bX, N_4D * sizeof(double));
		cudaMemcpy(bkg_d.bX, m_bX.get_data().data(), N_4D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.bY, N_4D * sizeof(double));
		cudaMemcpy(bkg_d.bY, m_bY.get_data().data(), N_4D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.bZ, N_4D * sizeof(double));
		cudaMemcpy(bkg_d.bZ, m_bZ.get_data().data(), N_4D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.bmag, N_4D * sizeof(double));
		cudaMemcpy(bkg_d.bmag, m_bmag.get_data().data(), N_4D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.eX, N_4D * sizeof(double));
		cudaMemcpy(bkg_d.eX, m_eX.get_data().data(), N_4D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.eY, N_4D * sizeof(double));
		cudaMemcpy(bkg_d.eY, m_eY.get_data().data(), N_4D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.eZ, N_4D * sizeof(double));
		cudaMemcpy(bkg_d.eZ, m_eZ.get_data().data(), N_4D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.emag, N_4D * sizeof(double));
		cudaMemcpy(bkg_d.emag, m_emag.get_data().data(), N_4D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.uX, N_4D * sizeof(double));
		cudaMemcpy(bkg_d.uX, m_uX.get_data().data(), N_4D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.uY, N_4D * sizeof(double));
		cudaMemcpy(bkg_d.uY, m_uY.get_data().data(), N_4D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.uZ, N_4D * sizeof(double));
		cudaMemcpy(bkg_d.uZ, m_uZ.get_data().data(), N_4D * sizeof(double), 
			cudaMemcpyHostToDevice);

		// 3D arrays
		cudaMalloc(&bkg_d.dXdx, N_3D * sizeof(double));
		cudaMemcpy(bkg_d.dXdx, m_dXdx.get_data().data(), N_3D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.dYdx, N_3D * sizeof(double));
		cudaMemcpy(bkg_d.dYdx, m_dYdx.get_data().data(), N_3D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.dZdx, N_3D * sizeof(double));
		cudaMemcpy(bkg_d.dZdx, m_dZdx.get_data().data(), N_3D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.dXdy, N_3D * sizeof(double));
		cudaMemcpy(bkg_d.dXdy, m_dXdy.get_data().data(), N_3D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.dYdy, N_3D * sizeof(double));
		cudaMemcpy(bkg_d.dYdy, m_dYdy.get_data().data(), N_3D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.dZdy, N_3D * sizeof(double));
		cudaMemcpy(bkg_d.dZdy, m_dZdy.get_data().data(), N_3D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.dXdz, N_3D * sizeof(double));
		cudaMemcpy(bkg_d.dXdz, m_dXdz.get_data().data(), N_3D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.dYdz, N_3D * sizeof(double));
		cudaMemcpy(bkg_d.dYdz, m_dYdz.get_data().data(), N_3D * sizeof(double), 
			cudaMemcpyHostToDevice);

		cudaMalloc(&bkg_d.dZdz, N_3D * sizeof(double));
		cudaMemcpy(bkg_d.dZdz, m_dZdz.get_data().data(), N_3D * sizeof(double), 
			cudaMemcpyHostToDevice);

		return bkg_d;

#else
		std::cerr << "Error! Background::to_device() was called but GPU support"
				  << " was not compiled in.\n";
#endif
	}  // to_device


	// Free up memory from device-side BackgroundDevice object
	void free_bkg(BackgroundDevice& bkg_d, int device_id)
	{

#ifdef USE_CUDA

		cudaSetDevice(device_id);

		// Free only the pointers that live in global memory. We don't free
		// the constant-memory arrays.

		cudaFree(bkg_d.ne);      bkg_d.ne = nullptr;
		cudaFree(bkg_d.te);      bkg_d.te = nullptr;
		cudaFree(bkg_d.ti);      bkg_d.ti = nullptr;
		//cudaFree(bkg_d.vp);      bkg_d.vp = nullptr;

		cudaFree(bkg_d.bX);      bkg_d.bX = nullptr;
		cudaFree(bkg_d.bY);      bkg_d.bY = nullptr;
		cudaFree(bkg_d.bZ);      bkg_d.bZ = nullptr;
		cudaFree(bkg_d.bmag);    bkg_d.bmag = nullptr;

		//cudaFree(bkg_d.gradbX); bkg_d.gradbX = nullptr;
		//cudaFree(bkg_d.gradbY); bkg_d.gradbY = nullptr;
		//cudaFree(bkg_d.gradbZ); bkg_d.gradbZ = nullptr;

		cudaFree(bkg_d.eX);      bkg_d.eX = nullptr;
		cudaFree(bkg_d.eY);      bkg_d.eY = nullptr;
		cudaFree(bkg_d.eZ);      bkg_d.eZ = nullptr;
		cudaFree(bkg_d.emag);    bkg_d.emag = nullptr;

		cudaFree(bkg_d.uX);      bkg_d.uX = nullptr;
		cudaFree(bkg_d.uY);      bkg_d.uY = nullptr;
		cudaFree(bkg_d.uZ);      bkg_d.uZ = nullptr;

		cudaFree(bkg_d.X);       bkg_d.X = nullptr;
		cudaFree(bkg_d.Y);       bkg_d.Y = nullptr;
		cudaFree(bkg_d.Z);       bkg_d.Z = nullptr;

		cudaFree(bkg_d.grid_X);  bkg_d.grid_X = nullptr;
		cudaFree(bkg_d.grid_Y);  bkg_d.grid_Y = nullptr;
		cudaFree(bkg_d.grid_Z);  bkg_d.grid_Z = nullptr;

		cudaFree(bkg_d.J);       bkg_d.J = nullptr;

		cudaFree(bkg_d.gij_00);  bkg_d.gij_00 = nullptr;
		cudaFree(bkg_d.gij_01);  bkg_d.gij_01 = nullptr;
		cudaFree(bkg_d.gij_02);  bkg_d.gij_02 = nullptr;
		cudaFree(bkg_d.gij_11);  bkg_d.gij_11 = nullptr;
		cudaFree(bkg_d.gij_12);  bkg_d.gij_12 = nullptr;
		cudaFree(bkg_d.gij_22);  bkg_d.gij_22 = nullptr;

		// dxdX/dydX/dzdX etc. are constant-memory, so we don't free

		cudaFree(bkg_d.dXdx);    bkg_d.dXdx = nullptr;
		cudaFree(bkg_d.dYdx);    bkg_d.dYdx = nullptr;
		cudaFree(bkg_d.dZdx);    bkg_d.dZdx = nullptr;

		cudaFree(bkg_d.dXdy);    bkg_d.dXdy = nullptr;
		cudaFree(bkg_d.dYdy);    bkg_d.dYdy = nullptr;
		cudaFree(bkg_d.dZdy);    bkg_d.dZdy = nullptr;

		cudaFree(bkg_d.dXdz);    bkg_d.dXdz = nullptr;
		cudaFree(bkg_d.dYdz);    bkg_d.dYdz = nullptr;
		cudaFree(bkg_d.dZdz);    bkg_d.dZdz = nullptr;

#endif

	}  // free_bkg
}
