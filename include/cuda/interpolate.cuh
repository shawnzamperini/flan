#pragma once

/**
* @file  interpolate.cuh
*
* @brief CUDA routines for interpolating background fields at a particle's
*        location.
*
* The background fields are defined at the grid cell centers (and, for 4D
* fields, at each time frame). To interpolate at a particle, a stencil is
* first built with build_stencil_3d or build_stencil_4d. It holds the
* indices of the 8 cell centers that bound the particle and the geometry
* of the box they form. One stencil can then be reused for every field
* interpolated at that particle with interpolate_field_3d or
* interpolate_field_4d, which saves recalculating the indices for each
* field. Everything is inline so it can be used in any kernel.
*
* The neighbors and centers come from the functions in zones_hd.h, the
* same ones the CPU uses, so the stencil is the CPU's (Boris::interp_recp
* and interp_tang).
*/
#include <cstddef>

#include "background_device.h"

#include "utilities.cuh"
#include "utilities_hd.h"
#include "zones_hd.h"


namespace Interpolate
{

	/**
	* Indices and geometry of the box of cell centers surrounding a point,
	* used to interpolate fields at that point.
	*
	* The same struct is used for 3D and 4D fields. A 3D stencil from
	* build_stencil_3d leaves the time members (base1, tidx_neighbor, t0
	* and frame_dt) unset.
	*
	* A corner's value is field[base0 + local[c]] (at tidx) or
	* field[base1 + local[c]] (at tidx_neighbor): base is where the zone's
	* frame starts in the flat ZonedVector data, and local is the corner's
	* index within it, which is small enough for an int.
	*
	* The 8 corners are ordered by (x, y, z) offset as 000, 001, 010, 011,
	* 100, 101, 110, 111, where 0 is the containing cell's index and 1 is
	* the neighbor's index in that direction.
	*/
	struct InterpolationStencil
	{
		std::size_t base0;   ///< Start of the zone's frame tidx (or of the
		                     ///< zone, for 3D).
		std::size_t base1;   ///< Start of the zone's frame tidx_neighbor
		                     ///< (4D only).
		int local[8];        ///< Corner indices within the zone's frame.

		int tidx_neighbor;   ///< Neighboring time frame index (4D only).
		int xidx_neighbor;   ///< Neighboring x cell index (may be a ghost).
		int yidx_neighbor;   ///< Neighboring y cell index (may be a ghost).
		int zidx_neighbor;   ///< Neighboring z cell index (may be a ghost).

		double t0;           ///< Time of frame tidx (s) (4D only).
		double frame_dt;     ///< Time from frame tidx to tidx_neighbor (s).
		                     ///< Negative if the neighbor is the previous
		                     ///< frame (4D only).

		double x0;           ///< x of the containing cell's center.
		double y0;           ///< y of the containing cell's center.
		double z0;           ///< z of the containing cell's center.

		double dx;           ///< x from the containing cell's center to the
		                     ///< neighbor's. Negative if the neighbor is
		                     ///< on the low side.
		double dy;           ///< Same as dx, for y.
		double dz;           ///< Same as dx, for z.
	};


	/**
	* Fills in the spatial part of a stencil: the neighbors, the corners'
	* local indices and the box of centers.
	*
	* @param zones The zones on the device.
	* @param zone Zone the point is in.
	* @param xidx x index of the grid cell containing the point.
	* @param yidx y index of the grid cell containing the point.
	* @param zidx z index of the grid cell containing the point.
	* @param x Computational x coordinate of the point.
	* @param y Computational y coordinate of the point.
	* @param z Computational z coordinate of the point.
	* @param[out] s Stencil to fill in.
	*/
	__device__ __forceinline__
	void build_stencil_space(const Zones::ZonesView& zones, const int zone,
		const int xidx, const int yidx, const int zidx, const double x,
		const double y, const double z, InterpolationStencil& s)
	{
		// The neighbor on the point's side of the cell center, which at the
		// edge of a zone may be a ghost cell across a shared face
		const int xn {Zones::neighbor_index(zones, zone, 0, x, xidx)};
		const int yn {Zones::neighbor_index(zones, zone, 1, y, yidx)};
		const int zn {Zones::neighbor_index(zones, zone, 2, z, zidx)};
		s.xidx_neighbor = xn;
		s.yidx_neighbor = yn;
		s.zidx_neighbor = zn;

		// Strides of the zone's frame, ghost cells included (as
		// Zones::flat_index)
		const int g {zones.ghost};
		const int pz {zones.n[2][zone] + 2 * g};
		const int pyz {(zones.n[1][zone] + 2 * g) * pz};
		const int i0 {(xidx + g) * pyz};
		const int i1 {(xn + g) * pyz};
		const int j0 {(yidx + g) * pz};
		const int j1 {(yn + g) * pz};
		const int k0 {zidx + g};
		const int k1 {zn + g};
		s.local[0] = i0 + j0 + k0;
		s.local[1] = i0 + j0 + k1;
		s.local[2] = i0 + j1 + k0;
		s.local[3] = i0 + j1 + k1;
		s.local[4] = i1 + j0 + k0;
		s.local[5] = i1 + j0 + k1;
		s.local[6] = i1 + j1 + k0;
		s.local[7] = i1 + j1 + k1;

		// The box between the cell centers (not the cell edges), since
		// that's where the fields are defined. A ghost cell's center is in
		// this zone's frame.
		s.x0 = Zones::center(zones, zone, 0, xidx);
		s.y0 = Zones::center(zones, zone, 1, yidx);
		s.z0 = Zones::center(zones, zone, 2, zidx);
		s.dx = Zones::center(zones, zone, 0, xn) - s.x0;
		s.dy = Zones::center(zones, zone, 1, yn) - s.y0;
		s.dz = Zones::center(zones, zone, 2, zn) - s.z0;
	}


	/**
	* Builds the stencil used to interpolate a 3D (x, y, z) field, such as
	* the basis vectors, at a point.
	*
	* @param bkg_d Background plasma on the device.
	* @param zone Zone the point is in.
	* @param xidx x index of the grid cell containing the point.
	* @param yidx y index of the grid cell containing the point.
	* @param zidx z index of the grid cell containing the point.
	* @param x Computational x coordinate of the point.
	* @param y Computational y coordinate of the point.
	* @param z Computational z coordinate of the point.
	* @return Stencil for interpolate_field_3d. The time members are unset.
	*/
	__device__ __forceinline__
	InterpolationStencil build_stencil_3d(
		const Background::BackgroundDevice& bkg_d, const int zone, int xidx,
		int yidx, int zidx, const double x, const double y, const double z)
	{
		InterpolationStencil s;
		build_stencil_space(bkg_d.zones, zone, xidx, yidx, zidx, x, y, z, s);
		s.base0 = bkg_d.zones.cell_offset[zone];
		return s;
	}


	/**
	* Builds the stencil used to interpolate a 4D (t, x, y, z) field, such
	* as E, B, ne or Te, at a point.
	*
	* The neighbors in space are picked as in build_stencil_3d, and in time
	* the neighboring frame is the one on the same side of frame tidx as t
	* (or the only neighbor at the first or last frame), so t is bounded by
	* the two frames.
	*
	* @param bkg_d Background plasma on the device.
	* @param zone Zone the point is in.
	* @param tidx Index of the nearest background time frame.
	* @param xidx x index of the grid cell containing the point.
	* @param yidx y index of the grid cell containing the point.
	* @param zidx z index of the grid cell containing the point.
	* @param t Time of the point (s).
	* @param x Computational x coordinate of the point.
	* @param y Computational y coordinate of the point.
	* @param z Computational z coordinate of the point.
	* @return Stencil for interpolate_field_4d.
	*/
	__device__ __forceinline__
	InterpolationStencil build_stencil_4d(
		const Background::BackgroundDevice& bkg_d, const int zone, int tidx,
		int xidx, int yidx, int zidx, const double t, const double x,
		const double y, const double z)
	{
		InterpolationStencil s;
		build_stencil_space(bkg_d.zones, zone, xidx, yidx, zidx, x, y, z, s);

		const int tn {Utilities::get_neighbor_index_cuda(t, bkg_d.times,
			tidx, bkg_d.tdim)};
		s.tidx_neighbor = tn;

		// Zone z's frame t starts at tdim * cell_offset[z] + t * (cells in
		// the zone), as in Zones::flat_index
		const std::size_t zone_start {static_cast<std::size_t>(bkg_d.tdim)
			* bkg_d.zones.cell_offset[zone]};
		const std::size_t frame {bkg_d.zones.cell_offset[zone + 1]
			- bkg_d.zones.cell_offset[zone]};
		s.base0 = zone_start + static_cast<std::size_t>(tidx) * frame;
		s.base1 = zone_start + static_cast<std::size_t>(tn) * frame;

		s.t0 = bkg_d.times[tidx];
		s.frame_dt = bkg_d.times[tn] - s.t0;
		return s;
	}


	/**
	* Trilinearly interpolates a 3D field at (x, y, z).
	*
	* The field values at the 8 stencil corners are passed to
	* trilinear_box (utilities_hd.h, shared with the CPU), which clamps
	* the weights to [0, 1], so there
	* is no extrapolation past the box of cell centers.
	*
	* @param field Flattened 3D field on the device.
	* @param s Stencil from build_stencil_3d, built at this point.
	* @param x Computational x coordinate of the point.
	* @param y Computational y coordinate of the point.
	* @param z Computational z coordinate of the point.
	* @return Interpolated field value.
	*/
	__device__ __forceinline__
	double interpolate_field_3d(const double* __restrict__ field,
		const InterpolationStencil& s, double x, double y, double z)
	{

		// Values at the 8 corners of the spatial interpolation stencil
		double verts0[8];

		// Calculate the 8 indices at tidx. 
		// pragma unroll just tells the compiler to break this loop 
		// down into 8 individual statements. This lets it use the 
		// register more efficiently and removes loop overhead.
		#pragma unroll
		for (int i = 0; i < 8; ++i)
		{
			verts0[i] = field[s.base0 + s.local[i]];
		}

		// Trilinear interpolation within the spatial cell defined
		// by the stencil.
		return  Utilities::trilinear_box(s.x0, s.y0, s.z0,
				s.dx, s.dy, s.dz, verts0[0], verts0[4], verts0[2], verts0[6],
				verts0[1], verts0[5], verts0[3], verts0[7], x, y, z);

	} // interpolate_field_3d


	/**
	* Interpolates a 4D field at (t, x, y, z), trilinear in space and
	* linear in time.
	*
	* The field is trilinearly interpolated at frames tidx and
	* tidx_neighbor as in interpolate_field_3d, and the two values are then
	* linearly interpolated to t. Unlike the spatial weights, the time
	* weight (t - t0) / frame_dt is not clamped. It is in [0, 1] except
	* before the first frame or after the last one, where the result is
	* extrapolated.
	*
	* @param field Flattened 4D field on the device.
	* @param s Stencil from build_stencil_4d, built at this point.
	* @param t Time of the point (s).
	* @param x Computational x coordinate of the point.
	* @param y Computational y coordinate of the point.
	* @param z Computational z coordinate of the point.
	* @return Interpolated field value.
	*/
	__device__ __forceinline__
	double interpolate_field_4d(const double* __restrict__ field,
		const InterpolationStencil& s, double t, double x, double y, double z)
	{

		// Arrays to hold all 8 values at the vertices
		double verts0[8];
		double verts1[8];

		// Calculate the 8 indices at tidx and tidx_neighbor. 
		// pragma unroll just tells the compiler to break this loop 
		// down into 8x2 individual statements. This lets it use the 
		// register more efficiently and removes loop overhead.
		#pragma unroll
		for (int i = 0; i < 8; ++i)
		{
			verts0[i] = field[s.base0 + s.local[i]];
			verts1[i] = field[s.base1 + s.local[i]];
		}

		// Trilinear interpolation at the particle location for this
		// field component, comp0, and then comp1 is the same just
		// at the next time index. I've lined the indices of v up with
		// the corresponding ones from the trilinear_box 
		// signature.
		const double comp0 = Utilities::trilinear_box(s.x0, s.y0, s.z0,
				s.dx, s.dy, s.dz, verts0[0], verts0[4], verts0[2], verts0[6],
				verts0[1], verts0[5], verts0[3], verts0[7], x, y, z);

		const double comp1 = Utilities::trilinear_box(s.x0, s.y0, s.z0,
				s.dx, s.dy, s.dz, verts1[0], verts1[4], verts1[2], verts1[6],
				verts1[1], verts1[5], verts1[3], verts1[7], x, y, z);

		// Now linearly interpolate in time (this is just point slope)
		return comp0 + (comp1 - comp0) * ((t - s.t0) / s.frame_dt);

	} // interpolate_field_4d
}
