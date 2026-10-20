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
*/
#include "background_device.h"

#include "device_constants.cuh"
#include "utilities.cuh"


namespace Interpolate
{

	/**
	* Indices and geometry of the box of cell centers surrounding a point,
	* used to interpolate fields at that point.
	*
	* The same struct is used for 3D and 4D fields. A 3D stencil from
	* build_stencil_3d leaves the time members (idx1, tidx_neighbor, t0
	* and frame_dt) unset.
	*
	* The 8 corners are ordered by (x, y, z) offset as 000, 001, 010, 011,
	* 100, 101, 110, 111, where 0 is the containing cell's index and 1 is
	* the neighbor's index in that direction.
	*/
	struct InterpolationStencil
	{
		int idx0[8];         ///< Flattened indices of the corners at tidx.
		int idx1[8];         ///< Flattened indices of the corners at
		                     ///< tidx_neighbor (4D only).

		int tidx_neighbor;   ///< Neighboring time frame index (4D only).
		int xidx_neighbor;   ///< Neighboring x cell index.
		int yidx_neighbor;   ///< Neighboring y cell index.
		int zidx_neighbor;   ///< Neighboring z cell index.

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
	* Builds the stencil used to interpolate a 3D (x, y, z) field, such as
	* the basis vectors, at a point.
	*
	* For each direction, get_neighbor_index_cuda picks the neighboring
	* cell on the same side of the containing cell's center as the point
	* (or the only neighbor at a grid edge). The 8 corner indices are then
	* flattened with calc_3d_index_cuda, and the box origin and size are
	* taken from the cell center coordinates in device_constants.cuh.
	*
	* @param bkg_d Background plasma on the device, used for the grid
	*        dimensions.
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
		const Background::BackgroundDevice& bkg_d, int xidx, 
		int yidx, int zidx, const double x, const double y, const double z)
	{
		InterpolationStencil s;

		// Neighbor indices
		const int xidx_neighbor = Utilities::get_neighbor_index_cuda(
				x, d_x, xidx, bkg_d.xdim);

		const int yidx_neighbor = Utilities::get_neighbor_index_cuda(
				y, d_y, yidx, bkg_d.ydim);

		const int zidx_neighbor = Utilities::get_neighbor_index_cuda(
				z, d_z, zidx, bkg_d.zdim);

		s.xidx_neighbor = xidx_neighbor;
		s.yidx_neighbor = yidx_neighbor;
		s.zidx_neighbor = zidx_neighbor;

		// Indices at each vertex
		s.idx0[0] = Utilities::calc_3d_index_cuda(bkg_d.ydim, bkg_d.zdim, 
			xidx, yidx, zidx);

		s.idx0[1] = Utilities::calc_3d_index_cuda(bkg_d.ydim, bkg_d.zdim,
			xidx, yidx, zidx_neighbor);

		s.idx0[2] = Utilities::calc_3d_index_cuda(bkg_d.ydim, bkg_d.zdim,
			xidx, yidx_neighbor, zidx);

		s.idx0[3] = Utilities::calc_3d_index_cuda(bkg_d.ydim, bkg_d.zdim,
			xidx, yidx_neighbor, zidx_neighbor);

		s.idx0[4] = Utilities::calc_3d_index_cuda(bkg_d.ydim, bkg_d.zdim,
			xidx_neighbor, yidx, zidx);

		s.idx0[5] = Utilities::calc_3d_index_cuda(bkg_d.ydim, bkg_d.zdim,
			xidx_neighbor, yidx, zidx_neighbor);

		s.idx0[6] = Utilities::calc_3d_index_cuda(bkg_d.ydim, bkg_d.zdim,
			xidx_neighbor, yidx_neighbor, zidx);

		s.idx0[7] = Utilities::calc_3d_index_cuda(bkg_d.ydim, bkg_d.zdim,
			xidx_neighbor, yidx_neighbor, zidx_neighbor);


		// x, y, z coordinates of two bounding vertices to interpolate 
		// between. Note these are not grid vertices, but rather are formed
		// by cell center coordinates since that's where the field is assumed
		// to be defined. It is essentially a cell shifted by dx/2, dy/2
		// and dz/2 if that helps.
		// d_x, d_y, d_z are from device_constants.cuh
		s.x0 = d_x[xidx];
		s.y0 = d_y[yidx];
		s.z0 = d_z[zidx];

		const double x1 = d_x[xidx_neighbor];
		const double y1 = d_y[yidx_neighbor];
		const double z1 = d_z[zidx_neighbor];

		s.dx = x1 - s.x0;
		s.dy = y1 - s.y0;
		s.dz = z1 - s.z0;

		return s;
	}


	/**
	* Builds the stencil used to interpolate a 4D (t, x, y, z) field, such
	* as E, B, ne or Te, at a point.
	*
	* The neighbors are picked the same way as in build_stencil_3d,
	* including in time: the neighboring frame is the one on the same side
	* of frame tidx as t (or the only neighbor at the first or last frame),
	* so t is bounded by the two frames. The 8 corner indices are flattened
	* with calc_4d_index_cuda at both tidx (idx0) and tidx_neighbor (idx1).
	*
	* @param bkg_d Background plasma on the device, used for the grid
	*        dimensions.
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
		const Background::BackgroundDevice& bkg_d, int tidx, int xidx, 
		int yidx, int zidx, const double t, const double x, const double y, 
		const double z)
	{
		InterpolationStencil s;

		// Neighbor indices
		const int tidx_neighbor =
			Utilities::get_neighbor_index_cuda(t, d_t, tidx, bkg_d.tdim);

		const int xidx_neighbor =
			Utilities::get_neighbor_index_cuda(x, d_x, xidx, bkg_d.xdim);

		const int yidx_neighbor =
			Utilities::get_neighbor_index_cuda(y, d_y, yidx, bkg_d.ydim);

		const int zidx_neighbor =
			Utilities::get_neighbor_index_cuda(z, d_z, zidx, bkg_d.zdim);

		s.tidx_neighbor = tidx_neighbor;
		s.xidx_neighbor = xidx_neighbor;
		s.yidx_neighbor = yidx_neighbor;
		s.zidx_neighbor = zidx_neighbor;

		// Indices at t = tidx
		s.idx0[0] = Utilities::calc_4d_index_cuda(
			bkg_d.xdim, bkg_d.ydim, bkg_d.zdim,
			tidx, xidx, yidx, zidx);

		s.idx0[1] = Utilities::calc_4d_index_cuda(
			bkg_d.xdim, bkg_d.ydim, bkg_d.zdim,
			tidx, xidx, yidx, zidx_neighbor);

		s.idx0[2] = Utilities::calc_4d_index_cuda(
			bkg_d.xdim, bkg_d.ydim, bkg_d.zdim,
			tidx, xidx, yidx_neighbor, zidx);

		s.idx0[3] = Utilities::calc_4d_index_cuda(
			bkg_d.xdim, bkg_d.ydim, bkg_d.zdim,
			tidx, xidx, yidx_neighbor, zidx_neighbor);

		s.idx0[4] = Utilities::calc_4d_index_cuda(
			bkg_d.xdim, bkg_d.ydim, bkg_d.zdim,
			tidx, xidx_neighbor, yidx, zidx);

		s.idx0[5] = Utilities::calc_4d_index_cuda(
			bkg_d.xdim, bkg_d.ydim, bkg_d.zdim,
			tidx, xidx_neighbor, yidx, zidx_neighbor);

		s.idx0[6] = Utilities::calc_4d_index_cuda(
			bkg_d.xdim, bkg_d.ydim, bkg_d.zdim,
			tidx, xidx_neighbor, yidx_neighbor, zidx);

		s.idx0[7] = Utilities::calc_4d_index_cuda(
			bkg_d.xdim, bkg_d.ydim, bkg_d.zdim,
			tidx, xidx_neighbor, yidx_neighbor, zidx_neighbor);

		// Indices at t = tidx_neighbor
		s.idx1[0] = Utilities::calc_4d_index_cuda(
			bkg_d.xdim, bkg_d.ydim, bkg_d.zdim,
			tidx_neighbor, xidx, yidx, zidx);

		s.idx1[1] = Utilities::calc_4d_index_cuda(
			bkg_d.xdim, bkg_d.ydim, bkg_d.zdim,
			tidx_neighbor, xidx, yidx, zidx_neighbor);

		s.idx1[2] = Utilities::calc_4d_index_cuda(
			bkg_d.xdim, bkg_d.ydim, bkg_d.zdim,
			tidx_neighbor, xidx, yidx_neighbor, zidx);

		s.idx1[3] = Utilities::calc_4d_index_cuda(
			bkg_d.xdim, bkg_d.ydim, bkg_d.zdim,
			tidx_neighbor, xidx, yidx_neighbor, zidx_neighbor);

		s.idx1[4] = Utilities::calc_4d_index_cuda(
			bkg_d.xdim, bkg_d.ydim, bkg_d.zdim,
			tidx_neighbor, xidx_neighbor, yidx, zidx);

		s.idx1[5] = Utilities::calc_4d_index_cuda(
			bkg_d.xdim, bkg_d.ydim, bkg_d.zdim,
			tidx_neighbor, xidx_neighbor, yidx, zidx_neighbor);

		s.idx1[6] = Utilities::calc_4d_index_cuda(
			bkg_d.xdim, bkg_d.ydim, bkg_d.zdim,
			tidx_neighbor, xidx_neighbor, yidx_neighbor, zidx);

		s.idx1[7] = Utilities::calc_4d_index_cuda(
			bkg_d.xdim, bkg_d.ydim, bkg_d.zdim,
			tidx_neighbor, xidx_neighbor,
			yidx_neighbor, zidx_neighbor);

		// t, x, y, z coordinates of two bounding vertices to interpolate 
		// between. Note these are not grid vertices, but rather are formed
		// by cell center coordinates since that's where the field is assumed
		// to be defined. It is essentially a cell shifted by dx/2, dy/2
		// and dz/2 if that helps.
		// d_t, d_x, d_y, d_z are from device_constants.cuh
		s.t0 = d_t[tidx];
		const double t1 = d_t[tidx_neighbor];

		s.x0 = d_x[xidx];
		s.y0 = d_y[yidx];
		s.z0 = d_z[zidx];

		const double x1 = d_x[xidx_neighbor];
		const double y1 = d_y[yidx_neighbor];
		const double z1 = d_z[zidx_neighbor];

		s.dx = x1 - s.x0;
		s.dy = y1 - s.y0;
		s.dz = z1 - s.z0;

		s.frame_dt = t1 - s.t0;

		return s;
	}


	/**
	* Trilinearly interpolates a 3D field at (x, y, z).
	*
	* The field values at the 8 stencil corners are passed to
	* trilinear_interpolate, which clamps the weights to [0, 1], so there
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
			verts0[i] = field[s.idx0[i]];
		}

		// Trilinear interpolation within the spatial cell defined
		// by the stencil.
		return  Utilities::trilinear_interpolate(s.x0, s.y0, s.z0,
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
			verts0[i] = field[s.idx0[i]];
			verts1[i] = field[s.idx1[i]];
		}

		// Trilinear interpolation at the particle location for this
		// field component, comp0, and then comp1 is the same just
		// at the next time index. I've lined the indices of v up with
		// the corresponding ones from the trilinear_interpolate 
		// signature.
		const double comp0 = Utilities::trilinear_interpolate(s.x0, s.y0, s.z0,
				s.dx, s.dy, s.dz, verts0[0], verts0[4], verts0[2], verts0[6],
				verts0[1], verts0[5], verts0[3], verts0[7], x, y, z);

		const double comp1 = Utilities::trilinear_interpolate(s.x0, s.y0, s.z0,
				s.dx, s.dy, s.dz, verts1[0], verts1[4], verts1[2], verts1[6],
				verts1[1], verts1[5], verts1[3], verts1[7], x, y, z);

		// Now linearly interpolate in time (this is just point slope)
		return comp0 + (comp1 - comp0) * ((t - s.t0) / s.frame_dt);

	} // interpolate_field_4d
}
