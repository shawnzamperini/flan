#pragma once
#include <cstdio>

#include "background_device.h"

#include "device_constants.cuh"
#include "interpolate.cuh"
#include "utilities.cuh"


namespace Interpolate
{
	
	// A stencil to hold the indices of cells surrounding a point that we
	// wish to interpolate at. We reuse this for 3D or 4D, just means in 3D
	// not all the fields are defined.
	struct InterpolationStencil
	{
		int idx0[8];     // corners at tidx
		int idx1[8];     // corners at tidx_neighbor

		int tidx_neighbor;
		int xidx_neighbor;
		int yidx_neighbor;
		int zidx_neighbor;

		double t0;
		double frame_dt;

		double x0;
		double y0;
		double z0;

		double dx;
		double dy;
		double dz;
	};


	// Create and return a stencil that is used in interpolate_field to 
	// interpolate a 3D field
	__device__ __forceinline__
	InterpolationStencil build_stencil_3d(
		const Background::BackgroundDevice& bkg_d, int xidx, 
		int yidx, int zidx)
	{
		InterpolationStencil s;

		// Neighbor indices
		const int xidx_neighbor = Utilities::get_neighbor_index_cuda(
				d_x[xidx], d_x, xidx, bkg_d.xdim);

		const int yidx_neighbor = Utilities::get_neighbor_index_cuda(
				d_y[yidx], d_y, yidx, bkg_d.ydim);

		const int zidx_neighbor = Utilities::get_neighbor_index_cuda(
				d_z[zidx], d_z, zidx, bkg_d.zdim);

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


	// Create and return a stencil that is used in interpolate_field to 
	// interpolate a 4D field
	__device__ __forceinline__
	InterpolationStencil build_stencil_4d(
		const Background::BackgroundDevice& bkg_d, int tidx, int xidx, 
		int yidx, int zidx)
	{
		InterpolationStencil s;

		// Neighbor indices
		const int xidx_neighbor =
			Utilities::get_neighbor_index_cuda(
				d_x[xidx], d_x, xidx, bkg_d.xdim);

		const int yidx_neighbor =
			Utilities::get_neighbor_index_cuda(
				d_y[yidx], d_y, yidx, bkg_d.ydim);

		const int zidx_neighbor =
			Utilities::get_neighbor_index_cuda(
				d_z[zidx], d_z, zidx, bkg_d.zdim);

		// Time neighbor
		const int at_end = (tidx == static_cast<int>(bkg_d.tdim) - 1);
		const int tidx_neighbor = tidx + 1 - 2 * at_end;

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


	// Interpolate a 3D field at x,y,z. The stencil comes from build_stencil
	// above. 
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


	// Interpolate a 4D field at t,x,y,z. The stencil comes from build_stencil
	// above. 
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
