/**
* @file  background.cu
*
* @brief GPU definitions of Background::to_device and free_bkg.
*
* The background is split between two kinds of GPU memory:
*
* - Constant memory: the times, cell centers, cell edges and simulation
*   bounds. These are small and are written to the __constant__ variables in
*   device_constants.cu. They don't belong to the BackgroundDevice.
* - Global memory: the plasma and field data and the basis vectors. These
*   are allocated with cudaMalloc, and their pointers are held in the
*   returned BackgroundDevice.
*
* CPU-only builds use the placeholder versions in background.cpp instead.
*/
#include <cuda_runtime.h>
#include <sstream>

#include "background.h"
#include "background_device.h"
#include "errors.h"

#include "device_constants.cuh"


namespace Background
{

	/**
	* Copies the background to a GPU.
	*
	* The times, cell centers, cell edges and simulation bounds are copied
	* to constant memory (see device_constants.cu), and the 3D and 4D arrays
	* are allocated and copied to global memory. Constant memory is per GPU,
	* so each GPU gets its own copy. Constant memory is shared by everything
	* on a GPU, though, so only one background can be on each GPU at a time.
	* Copying a second one to the same GPU overwrites the constant memory of
	* the first.
	*
	* Stops the run with Errors::fatal_error if a dimension is too big for
	* the constant memory arrays (see MAX_1D in device_constants.cuh).
	*
	* @param device_id ID of the GPU to copy to.
	* @return The device-side background. Its pointers are to global memory
	*         on the GPU and must be freed with free_bkg.
	*/
	BackgroundDevice Background::to_device(int device_id)
	{
		BackgroundDevice bkg_d {};

	#ifdef USE_CUDA

		bkg_d.device_id = device_id;

		// The times, cell centers and cell edges go in constant memory
		// arrays of size MAX_1D (see device_constants.cuh). If a dimension
		// doesn't fit, the copies below would write past the end of the
		// arrays, so stop the run. The cell edges have dim + 1 values, so
		// the spatial dimensions must be less than MAX_1D. dim1 is the time
		// dimension, which has no edges, so it can equal MAX_1D.
		if (m_dim1 > MAX_1D || m_dim2+1 > MAX_1D
			|| m_dim3+1 > MAX_1D || m_dim4+1 > MAX_1D)
		{
			std::ostringstream msg;
			msg << "MAX_1D is too small. Increase it in "
				<< "include/cuda/device_constants.cuh.\n"
				<< "  dim1 = " << m_dim1 << " (must be <= MAX_1D)\n"
				<< "  dim2 = " << m_dim2 << " (must be < MAX_1D)\n"
				<< "  dim3 = " << m_dim3 << " (must be < MAX_1D)\n"
				<< "  dim4 = " << m_dim4 << " (must be < MAX_1D)\n"
				<< "  MAX_1D = " << MAX_1D;
			Errors::fatal_error(msg.str());
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

		// These arrays are small and read by nearly every kernel, so they're
		// kept in constant memory, which is cached. They're defined in
		// device_constants.cu, and cudaMemcpyToSymbol copies into them on
		// the current GPU. The copies assume BkgFPType is double.
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

		// cudaMemcpyToSymbol copies from a host address, so for scalars we
		// pass the address of a local variable. The x, y and z bounds are
		// the outer cell edges, the same as the CPU uses.
		const double t_min {m_t_min};
		const double t_max {m_t_max};
		const double x_min {m_x_min};
		const double x_max {m_x_max};
		const double y_min {m_y_min};
		const double y_max {m_y_max};
		const double z_min {m_z_min};
		const double z_max {m_z_max};
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

		// Reciprocal basis vectors (grad x, grad y, grad z), 3D
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

		// Plasma and field data, 4D
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

		// Tangent basis vectors (dX/dx, etc.), 3D
		// These aren't actually used anywhere so commenting them out to save
		// memeory, but if they become useful just uncomment this section out.
		/*
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
		*/

		return bkg_d;

#else
		std::cerr << "Error! Background::to_device() was called but GPU support"
				  << " was not compiled in.\n";
#endif
	}  // to_device


	/**
	* Frees the global memory arrays of a device-side background.
	*
	* Each pointer is set to nullptr after it's freed. The constant memory
	* variables (d_t, d_x, etc.) are not freed, since constant memory is
	* allocated statically.
	*
	* @param bkg_d Device-side background from Background::to_device.
	* @param device_id ID of the GPU the background is on.
	*/
	void free_bkg(BackgroundDevice& bkg_d, int device_id)
	{

#ifdef USE_CUDA

		cudaSetDevice(device_id);

		// Free only the pointers that live in global memory. The constant
		// memory variables are static and can't be freed.

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

		cudaFree(bkg_d.dxdX);    bkg_d.dxdX = nullptr;
		cudaFree(bkg_d.dxdY);    bkg_d.dxdY = nullptr;
		cudaFree(bkg_d.dxdZ);    bkg_d.dxdZ = nullptr;

		cudaFree(bkg_d.dydX);    bkg_d.dydX = nullptr;
		cudaFree(bkg_d.dydY);    bkg_d.dydY = nullptr;
		cudaFree(bkg_d.dydZ);    bkg_d.dydZ = nullptr;

		cudaFree(bkg_d.dzdX);    bkg_d.dzdX = nullptr;
		cudaFree(bkg_d.dzdY);    bkg_d.dzdY = nullptr;
		cudaFree(bkg_d.dzdZ);    bkg_d.dzdZ = nullptr;

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
