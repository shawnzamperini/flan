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
*
* With the sheath on, its field (sheath_hd.h) is added to the background's,
* and particles that can reach it this step are pushed and moved in
* sub-steps (substep_kernel, the GPU version of Boris::substep_cpu).
*/
#include <cuda_runtime.h>
#include <cstdio>
#include <math.h>

#include "background_device.h"
#include "constants.h"
#include "sheath_hd.h"
#include "slots_device.h"
#include "zones_hd.h"

#include "interpolate.cuh"
#include "utilities.cuh"


namespace Boris
{

	/**
	* Advances the velocity of one particle by dt with the Boris algorithm,
	* the GPU version of push_one_cpu. The kernels (update_velocity_kernel,
	* substep_kernel) call it for their particles. The
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
	* In a zone with a zShift table, E (and |E|) are interpolated along the
	* field line instead (Interpolate::interp_plasma_vec3), as in
	* push_one_cpu.
	*
	* The sheath's field (Sheath::field_at) is added to E if sh's model
	* isn't off, as in push_one_cpu.
	*
	* @param slots_d Particle slots on the device. The Cartesian and
	*        computational velocities of particle i are updated in place.
	* @param bkg_d Background plasma on the device.
	* @param i Slot of the particle. Its cell indices must be those of its
	*        position.
	* @param dt Time step (s). Negative pushes backward.
	* @param sh The sheath on the device (model off for none).
	*/
	__device__ void push_one(const Slots::SlotsDevice& slots_d,
		const Background::BackgroundDevice& bkg_d, const int i,
		const double dt, const Sheath::SheathView& sh)
	{
		// Local variables
		int tidx {slots_d.tidx[i]};
		int xidx {slots_d.xidx[i]};
		int yidx {slots_d.yidx[i]};
		int zidx {slots_d.zidx[i]};
		const int zone {slots_d.zone[i]};

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
			{Interpolate::build_stencil_3d(bkg_d, zone, xidx, yidx, zidx, x, y, z)};
		Interpolate::InterpolationStencil stencil_4d
			{Interpolate::build_stencil_4d(bkg_d, zone, tidx, xidx, yidx, zidx, t, x, 
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

		// In a zone with a zShift table (a shifted-metric grid), E is
		// stretched along the field line like the rest of the plasma, so
		// it is interpolated along the field line instead, and rescaled to
		// the |E| interpolated the same way (as push_one_cpu). B is
		// axisymmetric and keeps the stencil's value.
		if (Zones::has_z_shift(bkg_d.zones, zone))
		{
			double e[3];
			Interpolate::interp_plasma_vec3(bkg_d, bkg_d.eX, bkg_d.eY,
				bkg_d.eZ, stencil_4d, zone, t, x, y, z, e);
			const double emag {Interpolate::interp_plasma(bkg_d, bkg_d.emag,
				stencil_4d, zone, t, x, y, z)};
			const double e_len {sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2])};
			const double scale {e_len > 0.0 ? emag / e_len : 0.0};
			for (int j {}; j < 3; ++j) E[j] = e[j] * scale;
			E[3] = emag;
		}

		// The sheath in front of a wall, far thinner than a cell, isn't
		// in the background's E: add its field in the cells next to
		// walls (E[3], the magnitude, isn't used by the push)
		if (sh.model != Sheath::Model::off)
		{
			const double pos[3] {x, y, z};
			double e_sh[3];
			Sheath::field_at(sh, zone, tidx, xidx, yidx, zidx, pos, e_sh);
			for (int j {}; j < 3; ++j) E[j] += e_sh[j];
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

	}  // push_one


	/**
	* Kernel that advances the velocity of each particle by one time step
	* with the Boris algorithm (push_one), one thread per slot. Dead
	* particles are skipped, and so are those sub-stepped near a wall this
	* step (nsub[i] > 1), which substep_kernel pushes instead.
	*
	* @param slots_d Particle slots on the device.
	* @param bkg_d Background plasma on the device.
	* @param dt Time step (s).
	* @param sh The sheath on the device (model off for none).
	* @param nsub Sub-steps of each slot, or nullptr for none.
	*/
	__global__ void update_velocity_kernel(Slots::SlotsDevice slots_d,
		const Background::BackgroundDevice bkg_d, const double dt,
		const Sheath::SheathView sh, const int* nsub)
	{
		const int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= slots_d.N) return;
		if (slots_d.state[i] > 0) return;
		if (nsub != nullptr && nsub[i] > 1) return;
		push_one(slots_d, bkg_d, i, dt, sh);
	}


	/**
	* Kernel that works out how many sub-steps each particle needs this
	* time step (Sheath::substeps): 1 for most, more for charged particles
	* that can reach the sheath field. As count_substeps_cpu.
	*
	* @param slots_d Particle slots on the device.
	* @param sh The sheath on the device.
	* @param dt Time step (s).
	* @param nsub Set to the sub-steps of each slot (1 for dead ones).
	*/
	__global__ void count_substeps_kernel(Slots::SlotsDevice slots_d,
		const Sheath::SheathView sh, const double dt, int* nsub)
	{
		const int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= slots_d.N) return;
		if (slots_d.state[i] > 0)
		{
			nsub[i] = 1;
			return;
		}
		const double q_m {slots_d.q[i] * -Constants::charge_e
			/ slots_d.mass};
		const double speed {sqrt(slots_d.vX[i] * slots_d.vX[i]
			+ slots_d.vY[i] * slots_d.vY[i] + slots_d.vZ[i]
			* slots_d.vZ[i])};
		const double pos[3] {slots_d.x[i], slots_d.y[i], slots_d.z[i]};
		nsub[i] = Sheath::substeps(sh, slots_d.zone[i], slots_d.tidx[i],
			slots_d.xidx[i], slots_d.yidx[i], slots_d.zidx[i], pos, speed,
			q_m, dt);
	}


	/**
	* Kernel that pushes the particles close to a wall through one time
	* step in sub-steps, one thread per slot: the GPU version of
	* Boris::substep_cpu (see it for the method). Particle i takes nsub[i]
	* sub-steps of h = dt / nsub[i] if nsub[i] > 1; the rest are left
	* alone. The positions are updated here; step_kernel still advances
	* the times.
	*
	* A thread with many sub-steps keeps its warp waiting; these are the
	* few particles within reach of the sheath.
	*
	* @param slots_d Particle slots on the device.
	* @param bkg_d Background plasma on the device.
	* @param dt The time step (s).
	* @param sh The sheath on the device.
	* @param nsub Sub-steps of each slot.
	*/
	__global__ void substep_kernel(Slots::SlotsDevice slots_d,
		const Background::BackgroundDevice bkg_d, const double dt,
		const Sheath::SheathView sh, const int* nsub)
	{
		const int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= slots_d.N) return;
		if (slots_d.state[i] > 0 || nsub[i] <= 1) return;

		const int n {nsub[i]};
		const double h {dt / n};
		const double t0 {slots_d.t[i]};
		const int zone {slots_d.zone[i]};
		const Zones::ZonesView& zones {bkg_d.zones};

		// The zone's outer cell edges, the grid a particle can leave
		double lo[3];
		double hi[3];
		for (int a {}; a < 3; ++a)
		{
			const double* g {Zones::zone_edges(zones, zone, a)};
			lo[a] = g[0];
			hi[a] = g[zones.n[a][zone]];
		}

		auto inside = [&]()
		{
			return slots_d.x[i] >= lo[0] && slots_d.x[i] <= hi[0]
				&& slots_d.y[i] >= lo[1] && slots_d.y[i] <= hi[1]
				&& slots_d.z[i] >= lo[2] && slots_d.z[i] <= hi[2];
		};
		auto set_cell = [&]()
		{
			slots_d.xidx[i] = Zones::cell_index(zones, zone, 0, slots_d.x[i]);
			slots_d.yidx[i] = Zones::cell_index(zones, zone, 1, slots_d.y[i]);
			slots_d.zidx[i] = Zones::cell_index(zones, zone, 2, slots_d.z[i]);
		};
		auto move = [&](const double tau)
		{
			slots_d.x[i] += slots_d.vx[i] * tau;
			slots_d.y[i] += slots_d.vy[i] * tau;
			slots_d.z[i] += slots_d.vz[i] * tau;
		};

		// v from t0 - dt/2 to t0 - h/2
		push_one(slots_d, bkg_d, i, 0.5 * (dt - h), sh);

		bool left {};
		for (int s {}; s < n; ++s)
		{
			slots_d.t[i] = t0 + s * h;
			if (s > 0) set_cell();
			push_one(slots_d, bkg_d, i, h, sh);
			const double before[3] {slots_d.x[i], slots_d.y[i],
				slots_d.z[i]};
			move(h);
			if (!inside())
			{
				// The fraction f of the sub-step at which it crossed the
				// grid's edge; a kick of (f - 1/2) h at the last position
				// inside gives its velocity when it crossed
				const double after[3] {slots_d.x[i], slots_d.y[i],
					slots_d.z[i]};
				double f {1.0};
				for (int a {}; a < 3; ++a)
				{
					const double step {after[a] - before[a]};
					if (after[a] < lo[a] && step < 0.0)
						f = fmin(f, (lo[a] - before[a]) / step);
					else if (after[a] > hi[a] && step > 0.0)
						f = fmin(f, (hi[a] - before[a]) / step);
				}
				f = fmin(fmax(f, 0.0), 1.0);
				slots_d.x[i] = before[0];
				slots_d.y[i] = before[1];
				slots_d.z[i] = before[2];
				slots_d.t[i] = t0 + (s + 0.5) * h;
				push_one(slots_d, bkg_d, i, (f - 0.5) * h, sh);
				slots_d.x[i] = after[0];
				slots_d.y[i] = after[1];
				slots_d.z[i] = after[2];
				move((n - s - 1) * h);
				left = true;
				break;
			}
		}

		// v from t0 + dt - h/2 back to t0 + dt/2
		if (!left)
		{
			slots_d.t[i] = t0 + dt;
			set_cell();
			push_one(slots_d, bkg_d, i, -0.5 * (dt - h), sh);
		}

		// step_kernel advances the time
		slots_d.t[i] = t0;
	}  // substep_kernel


	/**
	* Launches update_velocity_kernel to advance the velocity of every
	* particle on the GPU by one time step.
	*
	* Uses blocks of 256 threads, one thread per slot. The launch is
	* asynchronous, except in DEBUG builds, where it waits for the kernel to
	* finish and prints any error.
	*
	* With a sheath (sheath_d not nullptr) its field is added, and if
	* substep is true, the particles that can reach it are sub-stepped as
	* boris_wrapper does on the CPU: count_substeps_kernel sets
	* slots_d.nsub, update_velocity_kernel pushes the rest, and
	* substep_kernel pushes and moves these. step_gpu is then given
	* slots_d.nsub so it doesn't move them again. The half step back at
	* the start of the run adds the field but isn't sub-stepped, as on the
	* CPU.
	*
	* @param slots_d Particle slots on the device.
	* @param bkg_d Background plasma on the device.
	* @param dt Time step (s).
	* @param sheath_d The sheath on this GPU, or nullptr for none.
	* @param substep Sub-step the particles near the wall.
	*/
	void update_velocity_gpu(Slots::SlotsDevice& slots_d,
		const Background::BackgroundDevice& bkg_d, const double dt,
		const Sheath::SheathDevice* sheath_d, const bool substep)
	{
		
		// Block and grid size
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;

		const Sheath::SheathView sh {sheath_d != nullptr
			&& sheath_d->device_id >= 0 ? sheath_d->view
			: Sheath::SheathView {}};
		const bool sub {substep && sh.model != Sheath::Model::off};
		if (sub)
		{
			count_substeps_kernel<<<gridSize, blockSize>>>(slots_d, sh, dt,
				slots_d.nsub);
		}
		update_velocity_kernel<<<gridSize, blockSize>>>(slots_d, bkg_d, dt,
			sh, sub ? slots_d.nsub : nullptr);
		if (sub)
		{
			substep_kernel<<<gridSize, blockSize>>>(slots_d, bkg_d, dt, sh,
				slots_d.nsub);
		}

#ifdef DEBUG
		// Check for errors
		cudaError_t err {cudaDeviceSynchronize()};
		if (err != cudaSuccess)
			printf("update_velocity_kernel error: %s\n", cudaGetErrorString(err));
#endif

	} // update_velocity_gpu

} // namespace Boris
