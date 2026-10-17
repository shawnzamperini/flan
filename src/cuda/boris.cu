/**
* @file  boris.cu
*
* @brief CUDA routines for calculating the Lorentz force on a particle and
*        updating its velocity with the Boris algorithm.
*
* GPU version of boris.cpp. The velocity is updated in Cartesian space
* (vX, vY, vZ) using the electric and magnetic fields interpolated to the
* particle location. The velocity in computational coordinates (vx, vy, vz)
* is then found from the reciprocal basis vectors, which is what is used to
* move the particle.
*/
#include <cuda_runtime.h>
#include <cstdio>
#include <math.h>

#include "background_device.h"
#include "constants.h"
#include "slots_device.h"

#include "interpolate.cuh"
#include "utilities.cuh"


namespace Boris
{

	/**
	* Kernel that advances the velocity of each particle by one time step
	* using the Boris algorithm.
	*
	* One thread handles one slot, and dead particles are skipped. The
	* electric and magnetic fields (components and magnitudes) are
	* interpolated to the particle location, trilinear in space and linear
	* in time, with the 4D stencil from interpolate.cuh. This keeps the
	* fields continuous at cell boundaries and between time frames, which
	* avoids numerical diffusion and unphysical kicks. Each interpolated
	* vector is then rescaled to the interpolated magnitude so the
	* gyro-orbit is the correct size. A field whose interpolated vector is
	* zero is left as zero.
	*
	* The standard Boris rotation then updates the Cartesian velocity
	* (vX, vY, vZ). Velocities are stored at the half time steps before the
	* position, i.e., x_i = x(t_i) and v_i = v(t_i - dt/2). Finally, the
	* velocity in computational coordinates (vx, vy, vz) is calculated with
	* the reciprocal basis vectors, which are interpolated in space only
	* with the 3D stencil and are not rescaled.
	*
	* @param slots_d Particle slots on the device. The Cartesian and
	*        computational velocities are updated in place.
	* @param bkg_d Background plasma on the device.
	* @param dt Time step (s).
	*/
	__global__ void update_velocity_kernel(Slots::SlotsDevice slots_d, 
		const Background::BackgroundDevice bkg_d, const double dt)
	{
		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;

		// Don't try and access beyond the number of slots (segfault)
		if (i >= slots_d.N) return;

		// If particle isn't alive then skip
		if (slots_d.state[i] > 0) return;

		// Local variables
		int tidx {slots_d.tidx[i]};
		int xidx {slots_d.xidx[i]};
		int yidx {slots_d.yidx[i]};
		int zidx {slots_d.zidx[i]};

		double t {slots_d.t[i]};
		double x {slots_d.x[i]};
		double y {slots_d.y[i]};
		double z {slots_d.z[i]};

		// Trilinear interpolation to ensure we use a continuous B/E in the
		// algorithm. Without interpolating, the particle will:
		// a) Fail to follow the parallel direction of the field line due to
		//    it being discrete in the code (this is called "numerical
		//    diffusion" and introduces artificial cross-field transport).
		// b) Experience large "kicks" at cell boundaries where B/E is 
		//    discontinuous, which can introduce artificial drifts.
		// 
		// The issue in a) generally will still occur since we are still using
		// an approximation to the field line, but this greatly reduces the
		// numerical diffusion since B is now continuous.

		// Create stencils for interpolation
		Interpolate::InterpolationStencil stencil_3d
			{Interpolate::build_stencil_3d(bkg_d, xidx, yidx, zidx, x, y, z)};
		Interpolate::InterpolationStencil stencil_4d
			{Interpolate::build_stencil_4d(bkg_d, tidx, xidx, yidx, zidx, t, x, 
			y, z)};

		// Interpolate reciprocal basis vectors
		// 0-2: dxdX, dxdY, dxdZ
		// 3-5: dydX, dydY, dydZ
		// 6-8: dzdX, dzdY, dzdZ
		double rbv[9];
		rbv[0] = Interpolate::interpolate_field_3d(bkg_d.dxdX, stencil_3d, x, y, z);
		rbv[1] = Interpolate::interpolate_field_3d(bkg_d.dxdY, stencil_3d, x, y, z);
		rbv[2] = Interpolate::interpolate_field_3d(bkg_d.dxdZ, stencil_3d, x, y, z);
		rbv[3] = Interpolate::interpolate_field_3d(bkg_d.dydX, stencil_3d, x, y, z);
		rbv[4] = Interpolate::interpolate_field_3d(bkg_d.dydY, stencil_3d, x, y, z);
		rbv[5] = Interpolate::interpolate_field_3d(bkg_d.dydZ, stencil_3d, x, y, z);
		rbv[6] = Interpolate::interpolate_field_3d(bkg_d.dzdX, stencil_3d, x, y, z);
		rbv[7] = Interpolate::interpolate_field_3d(bkg_d.dzdY, stencil_3d, x, y, z);
		rbv[8] = Interpolate::interpolate_field_3d(bkg_d.dzdZ, stencil_3d, x, y, z);

		// Interpolate E and B components and magnitudes
		double E[4];
		E[0] = Interpolate::interpolate_field_4d(bkg_d.eX, stencil_4d, t, x, y, z);
		E[1] = Interpolate::interpolate_field_4d(bkg_d.eY, stencil_4d, t, x, y, z);
		E[2] = Interpolate::interpolate_field_4d(bkg_d.eZ, stencil_4d, t, x, y, z);
		E[3] = Interpolate::interpolate_field_4d(bkg_d.emag, stencil_4d, t, x, y, z);

		double B[4];
		B[0] = Interpolate::interpolate_field_4d(bkg_d.bX, stencil_4d, t, x, y, z);
		B[1] = Interpolate::interpolate_field_4d(bkg_d.bY, stencil_4d, t, x, y, z);
		B[2] = Interpolate::interpolate_field_4d(bkg_d.bZ, stencil_4d, t, x, y, z);
		B[3] = Interpolate::interpolate_field_4d(bkg_d.bmag, stencil_4d, t, x, y, z);

		// We interpolated the magnitude as well. This is because there is
		// no guarantee that the interpolated vector will have the same
		// magnitude (in fact it almost certainly shrinks). We are faced
		// with a choice:
		//   a) Leave as-is. B/E are continuous at cell boundaries, but
		//      the magnitude will be off. The size of the gyro-orbit will
		//      be affected as a result.
		//   b) Scale via the interpolated magnitude. B/E will be 
		//      discontinuous at cell boundaries, but definitely not as
		//      discontinuous as without interpolation. Gyro-orbit will
		//      be more physically correct.
		// We choose option b), because getting the gyro-orbit correct is
		// central to Flan. Slight discontinuities are not ideal, but the
		// particle spends order of magnitude more time not crossing 
		// a boundary, therefore option b) also minimizes numerical
		// diffusion in that sense as well by ensuring the gyro-orbit is
		// as correct as possible.
		//
		// We do not do this for the reciprocal basis vector. 
		// The reason is subtle and kind of tricky, but it has to do with
		// not destroying the orthogonality anymore than we have from
		// interpolating.

		// Compute magnitude of interpolated vectors
		double calc_Emag {sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2])};
		double calc_Bmag {sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2])};

		// Branch-free mask: 1.0 if calc_mag > 0, else 0.0
		double Emask = (calc_Emag > 0.0);
		double Bmask = (calc_Bmag > 0.0);

		// Scale vector components using interpolated magnitude. fmax prevents 
		// divide-by-zero
		double Edenom = fmax(calc_Emag, 1e-30);
		double Bdenom = fmax(calc_Bmag, 1e-30);
		for (int j {}; j < 3; ++j)
		{
			E[j] = E[j] * E[3] / Edenom * Emask;
			B[j] = B[j] * B[3] / Bdenom * Bmask;
		}

		// Store in local variables for convenience
		double q_m {slots_d.q[i] * -Constants::charge_e / slots_d.mass};
		double v[3] {slots_d.vX[i], slots_d.vY[i], slots_d.vZ[i]};

		// t vector (t is already taken so calling it tb, t-Boris, if you will)
		double tb[3];
		for (int j {}; j < 3; ++j) tb[j] = q_m * B[j] * 0.5 * dt;

		// Magnitude of t, squared
		double tmag2 {tb[0]*tb[0] + tb[1]*tb[1] + tb[2]*tb[2]};

		// s vector
		double s[3];
		for (int j {}; j < 3; ++j) s[j] = 2 * tb[j] / (1.0 + tmag2);

		// v minus
		double vminus[3];
		for (int j {}; j < 3; ++j) vminus[j] = v[j] + q_m * E[j] * 0.5 * dt;

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
		slots_d.vX[i] = vplus[0] + q_m * E[0] * 0.5 * dt;
		slots_d.vY[i] = vplus[1] + q_m * E[1] * 0.5 * dt;
		slots_d.vZ[i] = vplus[2] + q_m * E[2] * 0.5 * dt;

		// Update particle curvilinear velocity. For example,
		// vx = dx/dX * vX + dx/dY * vY + dx/dZ * vZ
		slots_d.vx[i] = rbv[0] * slots_d.vX[i] + rbv[1] * slots_d.vY[i] 
			+ rbv[2] * slots_d.vZ[i];
		slots_d.vy[i] = rbv[3] * slots_d.vX[i] + rbv[4] * slots_d.vY[i] 
			+ rbv[5] * slots_d.vZ[i];
		slots_d.vz[i] = rbv[6] * slots_d.vX[i] + rbv[7] * slots_d.vY[i] 
			+ rbv[8] * slots_d.vZ[i];

	}  // update_velocity_kernel


	/**
	* Launches update_velocity_kernel to advance the velocity of every
	* particle on the GPU by one time step.
	*
	* Uses blocks of 256 threads, one thread per slot. The launch is
	* asynchronous, except in DEBUG builds, where it waits for the kernel to
	* finish and prints any error.
	*
	* @param slots_d Particle slots on the device.
	* @param bkg_d Background plasma on the device.
	* @param dt Time step (s).
	*/
	void update_velocity_gpu(Slots::SlotsDevice& slots_d,
		const Background::BackgroundDevice& bkg_d, const double dt)
	{
		
		// Block and grid size
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;

		update_velocity_kernel<<<gridSize, blockSize>>>(slots_d, bkg_d, dt);

#ifdef DEBUG
		// Check for errors
		cudaError_t err {cudaDeviceSynchronize()};
		if (err != cudaSuccess)
			printf("update_velocity_kernel error: %s\n", cudaGetErrorString(err));
#endif

	} // update_velocity_gpu

} // namespace Boris
