/**
* @file  boundary.cu
*
* @brief CUDA routines for handling time and spatial boundary conditions.
*
* GPU version of check_bounds_cpu in boundary.cpp. Each time step, every
* particle is checked against the time boundary and the faces of its zone,
* with the per-particle routines shared with the CPU (boundary_hd.h): the
* time boundary, absorbing, wall, core, lcfs and periodic faces, crossings
* into other zones, the wall mask chi, and the turn of the velocity of a
* particle that goes around a toroidal wedge. What happened to the
* particles is counted on the device (add_fates_from_device reads it). The
* deposition on the wall and self-sputtering are CPU only so far.
*/
#include <cstdio>

#include "background_device.h"
#include "boundary.cuh"
#include "boundary_hd.h"
#include "errors.h"
#include "pcg32.h"
#include "slots_device.h"
#include "zones_hd.h"


namespace Boundary
{
	/**
	* Kernel that checks each particle against the time boundary and the
	* faces of its zone and applies the boundary conditions, one thread per
	* slot (see check_bounds_cpu for what happens at each kind of face).
	*
	* Dead particles are skipped. The time boundary comes first (time_bc),
	* then the faces (check_faces), then, if p.chi is set, the wall mask
	* (in_wall_cell), and last the turn of the Cartesian velocity of a
	* particle moved around the torus (turn_velocity). The counts of what
	* happened are added to fates_d with atomicAdd, only where nonzero,
	* since nearly every particle has nothing happen to it in a step.
	*
	* @param slots_d Particle slots on the device. The state, t, position,
	*        zone and (on a wedge) vX and vY are updated in place.
	* @param bkg_d Background plasma on the device, for the zones.
	* @param p Boundary settings.
	* @param rngs_d Random number generators on the device, one per slot,
	*        for the core faces.
	* @param fates_d Counters of what happened (num_fate_indices values),
	*        added to, or nullptr to not count.
	*/
	__global__
	void check_bounds_kernel(Slots::SlotsDevice slots_d,
		const Background::BackgroundDevice bkg_d, const BoundaryParams p,
		pcg32* rngs_d, unsigned long long* fates_d)
	{
		const int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= slots_d.N) return;
		if (slots_d.state[i] > 0) return;

		const Zones::ZonesView& zones {bkg_d.zones};
		FateTally tally {};

		// Time boundary
		double t {slots_d.t[i]};
		int state {time_bc(t, p, tally)};
		slots_d.t[i] = t;

		if (!state)
		{
			// Zone faces, on local copies
			double pos[3] {slots_d.x[i], slots_d.y[i], slots_d.z[i]};
			int zone {slots_d.zone[i]};
			double turn {};
			int wall_face {};
			state = check_faces(zones, p, pos, zone, turn, rngs_d[i], tally,
				wall_face);

			// Absorb particles that ended up in a cell inside the wall
			int idx[3];
			if (!state && p.chi != nullptr
				&& in_wall_cell(zones, p.chi, zone, pos, idx))
			{
				state = 1;
				++tally.n[fate_chi];
			}

			// The velocity is Cartesian, so it turns with the particle
			// around the torus
			if (turn != 0.0)
			{
				double vX {slots_d.vX[i]};
				double vY {slots_d.vY[i]};
				turn_velocity(vX, vY, p.phi_sign, turn);
				slots_d.vX[i] = vX;
				slots_d.vY[i] = vY;
			}

			slots_d.x[i] = pos[0];
			slots_d.y[i] = pos[1];
			slots_d.z[i] = pos[2];
			slots_d.zone[i] = zone;
		}
		slots_d.state[i] = state;

		if (fates_d != nullptr)
		{
			for (int k {}; k < num_fate_indices; ++k)
			{
				if (tally.n[k] != 0)
				{
					atomicAdd(&fates_d[k],
						static_cast<unsigned long long>(tally.n[k]));
				}
			}
		}
	}  // check_bounds_kernel


	/**
	* Launches check_bounds_kernel to apply the time and spatial boundary
	* conditions to every particle on the GPU.
	*
	* Uses blocks of 256 threads, one thread per slot. The launch is
	* asynchronous, except in DEBUG builds, where it waits for the kernel to
	* finish and prints any error.
	*
	* @param slots_d Particle slots on the device.
	* @param bkg_d Background plasma on the device.
	* @param p Boundary settings (see make_params). p.chi must be on the
	*        same GPU, e.g. bkg_d.chi.
	* @param rngs_d Random number generators on the device, one per slot.
	* @param fates_d Counters of what happened on this GPU
	*        (num_fate_indices values, from
	*        OpenADAS::IonizRecombWarningsDevice::fates), or nullptr.
	*/
	void check_bounds_gpu(Slots::SlotsDevice& slots_d,
		const Background::BackgroundDevice& bkg_d, const BoundaryParams& p,
		pcg32* rngs_d, unsigned long long* fates_d)
	{
		const int blockSize {256};
		const int gridSize {(slots_d.N + blockSize - 1) / blockSize};
		check_bounds_kernel<<<gridSize, blockSize>>>(slots_d, bkg_d, p,
			rngs_d, fates_d);

#ifdef DEBUG
		cudaError_t err {cudaDeviceSynchronize()};
		if (err != cudaSuccess)
			printf("check_bounds_kernel error: %s\n", cudaGetErrorString(err));
#endif
	}


	/**
	* Copies the counts of what happened to the particles on a GPU (from
	* check_bounds_kernel) and adds them to fates.
	*
	* Blocks until the GPU's kernels have finished.
	*
	* @param fates_d Counters on the device (num_fate_indices values).
	* @param device_id ID of the GPU they're on.
	* @param[in,out] fates Counts added to.
	*/
	void add_fates_from_device(const unsigned long long* fates_d,
		const int device_id, Fates& fates)
	{
		if (fates_d == nullptr) return;
		unsigned long long n[num_fate_indices] {};
		Errors::check_cuda(cudaSetDevice(device_id),
			"Boundary::add_fates_from_device");
		Errors::check_cuda(cudaMemcpy(n, fates_d, sizeof(n),
			cudaMemcpyDeviceToHost), "Boundary::add_fates_from_device");
		fates.crossings += static_cast<long long>(n[fate_crossings]);
		fates.absorbed += static_cast<long long>(n[fate_absorbed]);
		fates.wall += static_cast<long long>(n[fate_wall]);
		fates.chi += static_cast<long long>(n[fate_chi]);
		fates.core += static_cast<long long>(n[fate_core]);
		fates.time_end += static_cast<long long>(n[fate_time_end]);
		fates.runaway += static_cast<long long>(n[fate_runaway]);
	}

}  // namespace Boundary
