/**
* @file background_device.h
* @brief Contains struct to hold Background data for device-side access.
*/
#pragma once

#include "zones_hd.h"

namespace Background
{
	/**
	* Device-side copy of a Background, made by Background::to_device and
	* freed with free_bkg.
	*
	* The pointers are to global memory on the GPU. The fields are the flat
	* data of the host's ZonedVectors (zone after zone, ghost cells
	* included), so a value is found with Zones::flat_index on zones, with
	* dim1 = tdim for the 4D (t, x, y, z) fields and dim1 = 1 for the 3D
	* (x, y, z) ones. Vector components are in physical (Cartesian) X, Y and
	* Z.
	*/
	struct BackgroundDevice
	{
		int device_id;	///< ID of the GPU the arrays are on
		int tdim;		///< Number of background times
		double t_min;	///< First background time (s)
		double t_max;	///< Last background time (s)
		double* times;	///< Background times (s), tdim values

		/// The zones' grids, on the device (Zones::to_device)
		Zones::ZonesView zones;

		/// +1 if (R, z, Z) is right handed, -1 if not, for rotating
		/// vectors in toroidal zones (Background::get_phi_sign)
		double phi_sign;

		// Plasma and field data, 4D
		double* ne;		///< Electron density (m^-3)
		double* te;		///< Electron temperature (eV)
		double* ti;		///< Ion temperature (eV)
		//double* vp;

		double* bX;		///< Magnetic field X component (T)
		double* bY;		///< Magnetic field Y component (T)
		double* bZ;		///< Magnetic field Z component (T)
		double* bmag;	///< Magnetic field magnitude (T)

		//double* gradbX;
		//double* gradbY;
		//double* gradbZ;

		double* eX;		///< Electric field X component (V/m)
		double* eY;		///< Electric field Y component (V/m)
		double* eZ;		///< Electric field Z component (V/m)
		double* emag;	///< Electric field magnitude (V/m)

		double* uX;		///< Background flow X component (m/s)
		double* uY;		///< Background flow Y component (m/s)
		double* uZ;		///< Background flow Z component (m/s)

		/// Parallel ion and electron heat fluxes (W/m^2), for the thermal
		/// force, or nullptr if the background has none
		double* q_par_i;
		double* q_par_e;

		/// SOLEDGE3X wall mask chi, 3D (1 inside the wall), or nullptr if
		/// the background has none
		double* chi;

		// Not allocated by to_device, so these are always nullptr
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

		// Reciprocal (contravariant) basis vectors, 3D
		double* dxdX;	///< Reciprocal basis vector grad x, X component
		double* dxdY;	///< Reciprocal basis vector grad x, Y component
		double* dxdZ;	///< Reciprocal basis vector grad x, Z component

		double* dydX;	///< Reciprocal basis vector grad y, X component
		double* dydY;	///< Reciprocal basis vector grad y, Y component
		double* dydZ;	///< Reciprocal basis vector grad y, Z component

		double* dzdX;	///< Reciprocal basis vector grad z, X component
		double* dzdY;	///< Reciprocal basis vector grad z, Y component
		double* dzdZ;	///< Reciprocal basis vector grad z, Z component

		// Tangent (covariant) basis vectors, 3D
		double* dXdx;	///< Tangent basis vector dR/dx, X component
		double* dYdx;	///< Tangent basis vector dR/dx, Y component
		double* dZdx;	///< Tangent basis vector dR/dx, Z component

		double* dXdy;	///< Tangent basis vector dR/dy, X component
		double* dYdy;	///< Tangent basis vector dR/dy, Y component
		double* dZdy;	///< Tangent basis vector dR/dy, Z component

		double* dXdz;	///< Tangent basis vector dR/dz, X component
		double* dYdz;	///< Tangent basis vector dR/dz, Y component
		double* dZdz;	///< Tangent basis vector dR/dz, Z component

	};  // BackgroundDevice

}  // namespace Background
