/**
* @file  background.cu
*
* @brief GPU definitions of Background::to_device and free_bkg.
*
* The times, the zones' grids (Zones::to_device) and the plasma and field
* data and basis vectors are all copied to global memory, and their
* pointers are held in the returned BackgroundDevice. The fields are
* copied as the flat data of the host's ZonedVectors, ghost cells
* included.
*
* CPU-only builds use the placeholder versions in background.cpp instead.
*/
#include <cuda_runtime.h>
#include <string>

#include "background.h"
#include "background_device.h"
#include "errors.h"


namespace Background
{

	/**
	* Copies the background to a GPU.
	*
	* Everything goes to global memory: the times, the zones' grids (see
	* Zones::to_device) and the 3D and 4D arrays, as the flat data of their
	* ZonedVectors. Each GPU gets its own copy, and any number of
	* backgrounds can be on one GPU.
	*
	* Any background can be copied, but the GPU's push, collisions and
	* boundaries don't handle more than one zone or a toroidal wedge yet
	* (GPU plan G2 and G3), so a run stops before that (flan.cpp).
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

		Errors::check_cuda(cudaSetDevice(device_id), "Background::to_device");

		bkg_d.tdim = m_dim1;
		bkg_d.t_min = m_t_min;
		bkg_d.t_max = m_t_max;
		Errors::check_cuda(cudaMalloc(&bkg_d.times, m_dim1 * sizeof(double)),
			"Background::to_device: times");
		cudaMemcpy(bkg_d.times, m_times.data(), m_dim1 * sizeof(double),
			cudaMemcpyHostToDevice);
		bkg_d.zones = m_zones.to_device(device_id);

		// Sizes of the 3D and 4D arrays, ghost cells included (the last
		// zone offset is the total)
		const std::size_t N_3D {m_zones.view().cell_offset[
			m_zones.get_nzones()]};
		const std::size_t N_4D {static_cast<std::size_t>(m_dim1) * N_3D};

		// -----------------------------------
		// Global memory allocation and copies
		// -----------------------------------

		// Allocates a device array and copies a field into it. Every
		// field to_device copies must be loaded and shaped like the zones
		// (geometry with 1 frame, plasma with tdim), since the kernels
		// index them all with the same layout.
		auto copy_field = [](double*& dst, const auto& field,
			const std::size_t n, const char* name)
		{
			if (field.get_data().size() != n)
			{
				Errors::fatal_error(std::string("Background::to_device: ")
					+ name + " has " + std::to_string(field.get_data().size())
					+ " values, but the zones and times need "
					+ std::to_string(n) + ".");
			}
			Errors::check_cuda(cudaMalloc(&dst, n * sizeof(double)), name);
			Errors::check_cuda(cudaMemcpy(dst, field.get_data().data(),
				n * sizeof(double), cudaMemcpyHostToDevice), name);
		};

		// Reciprocal basis vectors (grad x, grad y, grad z), 3D
		copy_field(bkg_d.dxdX, m_dxdX, N_3D, "dxdX");
		copy_field(bkg_d.dxdY, m_dxdY, N_3D, "dxdY");
		copy_field(bkg_d.dxdZ, m_dxdZ, N_3D, "dxdZ");
		copy_field(bkg_d.dydX, m_dydX, N_3D, "dydX");
		copy_field(bkg_d.dydY, m_dydY, N_3D, "dydY");
		copy_field(bkg_d.dydZ, m_dydZ, N_3D, "dydZ");
		copy_field(bkg_d.dzdX, m_dzdX, N_3D, "dzdX");
		copy_field(bkg_d.dzdY, m_dzdY, N_3D, "dzdY");
		copy_field(bkg_d.dzdZ, m_dzdZ, N_3D, "dzdZ");

		// Plasma and field data, 4D
		copy_field(bkg_d.ne, m_ne, N_4D, "ne");
		copy_field(bkg_d.te, m_te, N_4D, "te");
		copy_field(bkg_d.ti, m_ti, N_4D, "ti");
		copy_field(bkg_d.bX, m_bX, N_4D, "bX");
		copy_field(bkg_d.bY, m_bY, N_4D, "bY");
		copy_field(bkg_d.bZ, m_bZ, N_4D, "bZ");
		copy_field(bkg_d.bmag, m_bmag, N_4D, "bmag");
		copy_field(bkg_d.eX, m_eX, N_4D, "eX");
		copy_field(bkg_d.eY, m_eY, N_4D, "eY");
		copy_field(bkg_d.eZ, m_eZ, N_4D, "eZ");
		copy_field(bkg_d.emag, m_emag, N_4D, "emag");
		copy_field(bkg_d.uX, m_uX, N_4D, "uX");
		copy_field(bkg_d.uY, m_uY, N_4D, "uY");
		copy_field(bkg_d.uZ, m_uZ, N_4D, "uZ");

		// The tangent basis vectors (dX/dx, etc.) aren't used by any
		// kernel, so they're left on the host to save memory

		return bkg_d;

#else
		std::cerr << "Error! Background::to_device() was called but GPU support"
				  << " was not compiled in.\n";
#endif
	}  // to_device


	/**
	* Frees the global memory arrays of a device-side background.
	*
	* Each pointer is set to nullptr after it's freed.
	*
	* @param bkg_d Device-side background from Background::to_device.
	* @param device_id ID of the GPU the background is on.
	*/
	void free_bkg(BackgroundDevice& bkg_d, int device_id)
	{

#ifdef USE_CUDA

		cudaSetDevice(device_id);

		cudaFree(bkg_d.times);   bkg_d.times = nullptr;
		Zones::free_zones(bkg_d.zones, device_id);

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
