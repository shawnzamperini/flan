#include <cuda_runtime.h>
#include <cstdio>
#include <math.h>

#include "background_device.h"
#include "constants.h"
#include "slots_device.h"

#include "device_constants.cuh"
#include "utilities.cuh"


namespace Boris
{

	__global__ void update_velocity_kernel(Slots::SlotsDevice slots_d, 
		const Background::BackgroundDevice bkg_d, const double dt)
	{
		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;

		// Don't try and access beyond the number of slots (segfault)
		if (i >= slots_d.N) return;

		// Local variables
		int tidx {slots_d.tidx[i]};
		int xidx {slots_d.xidx[i]};
		int yidx {slots_d.yidx[i]};
		int zidx {slots_d.zidx[i]};

		double t {slots_d.t[i]};
		double x {slots_d.x[i]};
		double y {slots_d.y[i]};
		double z {slots_d.z[i]};
		/*
		printf("x ptr = %p\n", slots_d.x);
		printf("y ptr = %p\n", slots_d.y);
		printf("z ptr = %p\n", slots_d.z);
		printf("i=%d\n", i);
		printf("tidx=%d xidx=%d yidx=%d zidx=%d\n", tidx, xidx, yidx, zidx);
		printf("t=%f x=%f y=%f z=%f\n", t, x, y, z);
		*/

		// Get nearest neighbor indices for each direction. These tell us
		// which direction we should interpolate towards, i.e., which
		// rectangle made by the neighboring cell centers our particle
		// is bounded by.
		// d_x is a global constant defined in cuda/device_constants.h.
		const int xidx_neighbor {Utilities::get_neighbor_index_cuda(x, 
			d_x, xidx, bkg_d.xdim)};
		const int yidx_neighbor {Utilities::get_neighbor_index_cuda(y, 
			d_y, yidx, bkg_d.ydim)};
		const int zidx_neighbor {Utilities::get_neighbor_index_cuda(z, 
			d_z, zidx, bkg_d.zdim)};

		// Similarly for t, except we can't use get_neighbor_index since it
		// uses cell center coordinates, and t is defined at each frame (i.e.,
		// not between each frame, that would be nonintuitive). So we use a
		// little SIMD-friendly logic here to assign tidx_neighbor to tidx+1,
		// and tidx-1 if we're in the last time frame.
		// at_end = 1 if tidx == ntimes-1, else 0 
		// If at_end = 0 → i + 1 
		// If at_end = 1 → i - 1 
		const int at_end {(tidx == static_cast<int>(bkg_d.tdim) - 1)}; 
		const int tidx_neighbor {tidx + 1 - 2 * at_end};

		// Trilinear interpolation to ensure we use a continous B/E in the 
		// algorithm. Without interpolating, the particle will:
		// a) Fail to follow the parallel direction of the field line due to
		//    it being discrete in the code (this is called "numerical
		//    diffusion" and introduces artificial cross-field transport).
		// b) Experience large "kicks" at cell boundaries where B/E is 
		//    discontinuous, which can introduce artificial drifts.
		// 
		// The issue in a) generally will still occur since we are still using
		// an approximation to the field line, but this greatly reduces the
		// numerical diffusion since B is now continous.

		// x, y, z coordinates of two bounding vertices to interpolate 
		// between. Note these are not grid vertices, but rather are formed
		// by cell center coordinates since that's where B/E are assumed
		// to be defined. It is essentially a cell shifted by dx/2, dy/2
		// and dz/2 if that helps.
		const double t0 {d_t[tidx]};
		const double t1 {d_t[tidx_neighbor]};
		const double x0 {d_x[xidx]};
		const double x1 {d_x[xidx_neighbor]};
		const double y0 {d_y[yidx]};
		const double y1 {d_y[yidx_neighbor]};
		const double z0 {d_z[zidx]};
		const double z1 {d_z[zidx_neighbor]};

		// Distance between neighbors
		const double dx {x1 - x0};
		const double dy {y1 - y0};
		const double dz {z1 - z0};

		// Time between neighboring t frames
		const double frame_dt {t1 - t0};

		// Array to hold field component and reciprocal basis vector pointers
		// 0: BX, BY, BZ, |B|
		// 1: EX, EY, EZ, |E|
		// 2: dxdX, dxdY, dxdZ, undefined
		// 3: dydX, dydY, dydZ, undefined
		// 4: dzdX, dzdY, dzdZ, undefined
		// We use restrict since these arrays never overlap in memory with
		// another pointer, and allows the compiler to optimize the arrays.
		double* __restrict__ field_ptrs[5][4];  

		// 4D arrays
		field_ptrs[0][0] = bkg_d.bX;
		field_ptrs[0][1] = bkg_d.bY;
		field_ptrs[0][2] = bkg_d.bZ;
		field_ptrs[0][3] = bkg_d.bmag;
		field_ptrs[1][0] = bkg_d.eX;
		field_ptrs[1][1] = bkg_d.eY;
		field_ptrs[1][2] = bkg_d.eZ;
		field_ptrs[1][3] = bkg_d.emag;

		// 3D arrays
		field_ptrs[2][0] = bkg_d.dxdX;
		field_ptrs[2][1] = bkg_d.dxdY;
		field_ptrs[2][2] = bkg_d.dxdZ;
		field_ptrs[3][0] = bkg_d.dydX;
		field_ptrs[3][1] = bkg_d.dydY;
		field_ptrs[3][2] = bkg_d.dydZ;
		field_ptrs[4][0] = bkg_d.dzdX;
		field_ptrs[4][1] = bkg_d.dzdY;
		field_ptrs[4][2] = bkg_d.dzdZ;

		// This will hold the interpolated field values
		double intrp_field[5][4];

		// 4D -> 1D index for indexing background arrays in bkg_d. 
		// _0 is at tidx, and _1 is at tidx_neighbor
		int idx000_0_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx, xidx, yidx, zidx)};
		int idx001_0_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx, xidx, yidx, zidx_neighbor)};
		int idx010_0_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx, xidx, yidx_neighbor, zidx)};
		int idx011_0_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx, xidx, yidx_neighbor, zidx_neighbor)};
		int idx100_0_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx, xidx_neighbor, yidx, zidx)};
		int idx101_0_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx, xidx_neighbor, yidx, zidx_neighbor)};
		int idx110_0_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx, xidx_neighbor, yidx_neighbor, zidx)};
		int idx111_0_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx, xidx_neighbor, yidx_neighbor, zidx_neighbor)};

		int idx000_1_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx_neighbor, xidx, yidx, zidx)};
		int idx001_1_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx_neighbor, xidx, yidx, zidx_neighbor)};
		int idx010_1_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx_neighbor, xidx, yidx_neighbor, zidx)};
		int idx011_1_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx_neighbor, xidx, yidx_neighbor, zidx_neighbor)};
		int idx100_1_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx_neighbor, xidx_neighbor, yidx, zidx)};
		int idx101_1_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx_neighbor, xidx_neighbor, yidx, zidx_neighbor)};
		int idx110_1_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx_neighbor, xidx_neighbor, yidx_neighbor, zidx)};
		int idx111_1_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx_neighbor, xidx_neighbor, yidx_neighbor, 
			zidx_neighbor)};

		// Again for the 3D->1D indices
		int idx000_0_3d {Utilities::calc_3d_index_cuda(bkg_d.ydim, 
			bkg_d.zdim, xidx, yidx, zidx)};
		int idx001_0_3d {Utilities::calc_3d_index_cuda(bkg_d.ydim, 
			bkg_d.zdim, xidx, yidx, zidx_neighbor)};
		int idx010_0_3d {Utilities::calc_3d_index_cuda(bkg_d.ydim, 
			bkg_d.zdim, xidx, yidx_neighbor, zidx)};
		int idx011_0_3d {Utilities::calc_3d_index_cuda(bkg_d.ydim, 
			bkg_d.zdim, xidx, yidx_neighbor, zidx_neighbor)};
		int idx100_0_3d {Utilities::calc_3d_index_cuda(bkg_d.ydim, 
			bkg_d.zdim, xidx_neighbor, yidx, zidx)};
		int idx101_0_3d {Utilities::calc_3d_index_cuda(bkg_d.ydim, 
			bkg_d.zdim, xidx_neighbor, yidx, zidx_neighbor)};
		int idx110_0_3d {Utilities::calc_3d_index_cuda(bkg_d.ydim, 
			bkg_d.zdim, xidx_neighbor, yidx_neighbor, zidx)};
		int idx111_0_3d {Utilities::calc_3d_index_cuda(bkg_d.ydim, 
			bkg_d.zdim, xidx_neighbor, yidx_neighbor, zidx_neighbor)};

		int idx000_1_3d {Utilities::calc_3d_index_cuda(bkg_d.ydim, 
			bkg_d.zdim, xidx, yidx, zidx)};
		int idx001_1_3d {Utilities::calc_3d_index_cuda(bkg_d.ydim, 
			bkg_d.zdim, xidx, yidx, zidx_neighbor)};
		int idx010_1_3d {Utilities::calc_3d_index_cuda(bkg_d.ydim, 
			bkg_d.zdim, xidx, yidx_neighbor, zidx)};
		int idx011_1_3d {Utilities::calc_3d_index_cuda(bkg_d.ydim, 
			bkg_d.zdim, xidx, yidx_neighbor, zidx_neighbor)};
		int idx100_1_3d {Utilities::calc_3d_index_cuda(bkg_d.ydim, 
			bkg_d.zdim, xidx_neighbor, yidx, zidx)};
		int idx101_1_3d {Utilities::calc_3d_index_cuda(bkg_d.ydim, 
			bkg_d.zdim, xidx_neighbor, yidx, zidx_neighbor)};
		int idx110_1_3d {Utilities::calc_3d_index_cuda(bkg_d.ydim, 
			bkg_d.zdim, xidx_neighbor, yidx_neighbor, zidx)};
		int idx111_1_3d {Utilities::calc_3d_index_cuda(bkg_d.ydim, 
			bkg_d.zdim, xidx_neighbor, yidx_neighbor, 
			zidx_neighbor)};

		// Group into an array, can help reduce register pressure
		int idx0_4d[8] = {idx000_0_4d, idx001_0_4d, idx010_0_4d, idx011_0_4d,
						idx100_0_4d, idx101_0_4d, idx110_0_4d, idx111_0_4d};
		int idx1_4d[8] = {idx000_1_4d, idx001_1_4d, idx010_1_4d, idx011_1_4d,
						idx100_1_4d, idx101_1_4d, idx110_1_4d, idx111_1_4d};

		int idx0_3d[8] = {idx000_0_3d, idx001_0_3d, idx010_0_3d, idx011_0_3d,
						idx100_0_3d, idx101_0_3d, idx110_0_3d, idx111_0_3d};
		int idx1_3d[8] = {idx000_1_3d, idx001_1_3d, idx010_1_3d, idx011_1_3d,
						idx100_1_3d, idx101_1_3d, idx110_1_3d, idx111_1_3d};

		// Array to hold all 8 values at the vertices
		double verts0[8];
		double verts1[8];

		// j = 0 = B, j = 1 = E, j = 2 = dx/dX,Y,Z 
		// j = 3 = dy/dX,Y,Z, j = 4 = dz/dX,Y,Z
		for (int j = 0; j < 5; ++j)
		{
			for (int k = 0; k < 4; ++k)
			{
				// We skip the last index for the reciprocal basis vector
				// components since it's garbage (we do not interpolate the
				// magnitude of the basis vectors because that would break
				// its properties, e.g., e_i * e^i = delta^i_j.
				if (j > 1 && k == 3) continue;

				// Calculate the 8 indices at tidx and tidx_neighbor. 
				// pragma unroll just tells the compiler to break this loop 
				// down into 8x2 individual statements. This lets it use the 
				// register more efficiently and removes loop overhead.
				#pragma unroll
				for (int l = 0; l < 8; ++l)
				{
					// 3D arrays
					if (j > 1)
					{
						verts0[l] = field_ptrs[j][k][idx0_3d[l]];
						verts1[l] = field_ptrs[j][k][idx1_3d[l]];
					}

					// 4D arrays
					else
					{
						verts0[l] = field_ptrs[j][k][idx0_4d[l]];
						verts1[l] = field_ptrs[j][k][idx1_4d[l]];
					}
					//if (i == 0 && j == 2)
					//	printf("%d %d verts[%d] = %f\n", j, k, l, verts0[l]);
				}

				// Trilinear interpolation at the particle location for this
				// field component, comp0, and then comp1 is the same just
				// at the next time index. I've lined the indices of v up with
				// the corresponding ones from the trilinear_interpolate 
				// signature.
				double comp0 {Utilities::trilinear_interpolate(x0, y0, z0, 
					dx, dy, dz, verts0[0], verts0[4], verts0[2], verts0[6], 
					verts0[1], verts0[5], verts0[3], verts0[7], x, y, z)};

				double comp1 {Utilities::trilinear_interpolate(x0, y0, z0, 
					dx, dy, dz, verts1[0], verts1[4], verts1[2], verts1[6], 
					verts1[1], verts1[5], verts1[3], verts1[7], x, y, z)};

				//if (i == 0) printf("%d comp0=%f, comp1=%f\n", j, comp0, comp1);

				// Now linearly interpolate in time (this is just point slope)
				const double slope {(comp1 - comp0) / frame_dt};
				intrp_field[j][k] = slope * (t - t0) + comp0;

			}  // k loop


			// We interpolated the magnitude as well. This is because there is
			// no guarantee that the interpolated vector will have the same
			// magnitude (in fact it almost certainly shrinks). We are faced
			// with a choice:
			//   a) Leave as-is. B/E are continuous at cell boundaries, but
			//      the magntiude will be off. The size of the gyro-orbit will
			//      be affected as a result.
			//   b) Scale via the interpolated magnitude. B/E will be 
			//      discontinuous at cell boundaries, but definitely not as
			//      discontinous as without interpolation. Gyro-orbit will
			//      be more physically correct.
			// We choose option b), because getting the gyro-orbit correct is
			// central to Flan. Slight discontinuities are not ideal, but the
			// particle spends order of magnitude more time not crossing 
			// a boundary, therefore option b) also minimizes numerical
			// diffusion in that sense as well by ensuring the gyro-orbit is
			// as correct as possible.

			// We do not do this for the reciprocal basis vector. 
			// The reason is subtle and kind of tricky, but it has to do with
			// not destroying the orthogonality anymore than we have from
			// interpolating.
			if (j > 1) continue;

			// Compute magnitude of interpolated vector
			double calc_mag {sqrt(intrp_field[j][0] * intrp_field[j][0] +
				intrp_field[j][1] * intrp_field[j][1] +
				intrp_field[j][2] * intrp_field[j][2])};

			// Branch-free mask: 1.0 if calc_mag > 0, else 0.0
			double mask = (calc_mag > 0.0);

			// Scale vector components using interpolated magnitude 
			// in interp_field[i][3]. fmax prevents divide-by-zero
			double denom = fmax(calc_mag, 1e-30);
			for (int k = 0; k < 3; ++k)
			{
				intrp_field[j][k] = intrp_field[j][k] * intrp_field[j][3] 
					/ denom * mask;
			}

		}  // j loop

		// Boris continues

		// Store in local variables for conveinence
		double q_m {slots_d.q[i] * -Constants::charge_e / slots_d.mass};
		double v[3] {slots_d.vX[i], slots_d.vY[i], slots_d.vZ[i]};

		// t vector (t is already taken so calling it tb, t-Boris, if you will)
		// Reminder that intrp_field[0] is B.
		double tb[3];
		for (int j {}; j < 3; ++j) tb[j] = q_m * intrp_field[0][j] * 0.5 * dt;

		// Magnitude of t, squared
		double tmag2 {tb[0]*tb[0] + tb[1]*tb[1] + tb[2]*tb[2]};

		// s vector
		double s[3];
		for (int j {}; j < 3; ++j) s[j] = 2 * tb[j] / (1.0 + tmag2);

		// v minus
		// intrp_field[1] is E
		double vminus[3];
		for (int j {}; j < 3; ++j) vminus[j] = v[j] + q_m * intrp_field[1][j] 
			* 0.5 * dt;

		// v prime
		double vprime[3];
		double vminus_cross_t[3]; 
		Utilities::cross_product_cuda(vminus, tb, vminus_cross_t);
		for (int j {}; j < 3; ++j) vprime[j] = vminus[j] + vminus_cross_t[j];

		// v plus
		double vplus[3];
		double vprime_cross_s[3]; 
		Utilities::cross_product_cuda(vprime, s, vprime_cross_s);
		for (int j {}; j < 3; ++j) vplus[j] = vminus[j] + vprime_cross_s[j];

		// v n+1/2
		// Note we are storing particle velocity at the half time steps before
		// the position. I.e., xi = x(ti), vi = v(ti - dt/2). 
		// intrp_field[1] is E
		slots_d.vX[i] = vplus[0] + q_m * intrp_field[1][0] * 0.5 * dt;
		slots_d.vY[i] = vplus[1] + q_m * intrp_field[1][1] * 0.5 * dt;
		slots_d.vZ[i] = vplus[2] + q_m * intrp_field[1][2] * 0.5 * dt;

		/*
		for (int j {}; j < 5; ++j)
		{
			for (int k {}; k < 4; ++k)
			{
				printf("interp_field[%d][%d] = %f\n", j, k, intrp_field[j][k]);
			}
		}
		*/

		//printf("slots_d.vX[i] = %f\n", slots_d.vX[i]);
		//printf("slots_d.vY[i] = %f\n", slots_d.vY[i]);
		//printf("slots_d.vZ[i] = %f\n", slots_d.vZ[i]);

		// Update particle curvilinear velocity. For example,
		// vx = dx/dX * vX + dx/dY * vY + dx/dZ * vZ
		slots_d.vx[i] = intrp_field[2][0] * slots_d.vX[i] 
			+ intrp_field[2][1] * slots_d.vY[i] 
			+ intrp_field[2][2] * slots_d.vZ[i];
		slots_d.vy[i] = intrp_field[3][0] * slots_d.vX[i] 
			+ intrp_field[3][1] * slots_d.vY[i] 
			+ intrp_field[3][2] * slots_d.vZ[i];
		slots_d.vz[i] = intrp_field[4][0] * slots_d.vX[i] 
			+ intrp_field[4][1] * slots_d.vY[i] 
			+ intrp_field[4][2] * slots_d.vZ[i];
		
	}  // update_velocity_kernel


	void update_velocity_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d, const double dt)
	{

		// Block and grid size
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;

		update_velocity_kernel<<<gridSize, blockSize>>>(slots_d, bkg_d, dt);

	} // update_velocity_gpu

} // namespace Boris
