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
* particles is counted on the device (add_fates_from_device reads it).
* Where absorbed atoms land on the wall is tallied on the device with the
* routines shared with the CPU (sputtering_hd.h), and added to the host's
* tallies after the run (Sputtering::Deposition::add_from_device). A
* particle landing on an element can carry on as an atom it sputtered
* (Sputtering::rebirth, the CPU's routine), with its slot's generator.
*/
#include <cstdio>

#include "background_device.h"
#include "boundary.cuh"
#include "boundary_hd.h"
#include "errors.h"
#include "pcg32.h"
#include "slots_device.h"
#include "sputtering_hd.h"
#include "zones_hd.h"


namespace Boundary
{
	/// Adds to the device's deposition tallies (atomicAdd, since many
	/// threads add to the same element or bin)
	struct DeviceAdd
	{
		__device__ void operator()(double* p, const double value) const
		{
			atomicAdd(p, value);
		}
	};

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
	* @param dep Deposition tallies on the device, added to where atoms
	*        land on the wall, if deposit is true.
	* @param deposit Whether to tally the deposition.
	* @param self Self-sputtering on the device, used if self_on is true
	*        (with deposit): a particle landing on an element can carry on
	*        as an atom it sputtered.
	* @param self_on Whether self-sputtering is on.
	* @param avalanche Set to 1 if a chain reached sput_self_max_gen.
	*/
	__global__
	void check_bounds_kernel(Slots::SlotsDevice slots_d,
		const Background::BackgroundDevice bkg_d, const BoundaryParams p,
		pcg32* rngs_d, unsigned long long* fates_d,
		const Sputtering::DepositionView dep, const bool deposit,
		const Sputtering::SelfView self, const bool self_on, int* avalanche)
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

			// Its velocity, for the impact energy and angle tallies (before
			// any turn, as on the CPU)
			const double v[3] {slots_d.vX[i], slots_d.vY[i], slots_d.vZ[i]};

			// Self-sputtering off element e, which the particle just landed
			// on at (pos[0], pos[1], z_wall), as check_bounds_cpu's
			// self_sputter: if it carries on as an atom it sputtered, the
			// slot gets that atom and true is returned
			auto self_sputter = [&](const long e, const double z_wall)
			{
				if (!self_on || e < 0) return false;
				const double hit[3] {pos[0], pos[1], z_wall};
				const Sputtering::Rebirth r {Sputtering::rebirth(self,
					static_cast<std::size_t>(e), slots_d.tidx[i], hit, v,
					slots_d.mass, slots_d.q[i], slots_d.weight[i],
					slots_d.generation[i], rngs_d[i], DeviceAdd {})};
				if (r.avalanche) atomicExch(avalanche, 1);
				if (!r.alive) return false;

				pos[0] = r.x;
				pos[1] = r.y;
				pos[2] = r.z;
				zone = r.zone;
				turn = 0.0;		// the new velocity is at the new z
				slots_d.xidx[i] = r.i;
				slots_d.yidx[i] = r.j;
				slots_d.zidx[i] = r.k;
				slots_d.vX[i] = r.vX;
				slots_d.vY[i] = r.vY;
				slots_d.vZ[i] = r.vZ;
				slots_d.vx[i] = 0.0;	// set by the next Boris update
				slots_d.vy[i] = 0.0;
				slots_d.vz[i] = 0.0;
				slots_d.q[i] = 0;
				slots_d.weight[i] = r.weight;
				slots_d.pt[i] = 0.0;	// time since its own sputtering
				slots_d.prompt_by[i] = -1.0;
				slots_d.generation[i] += 1;
				++tally.n[fate_self_sputtered];
				return true;
			};

			// Landed on a wall face: where on it, and maybe a
			// self-sputtered atom carries on from there
			if (state && wall_face >= 0 && deposit)
			{
				const long e {Sputtering::add_at_wall_face(dep, zones, zone,
					wall_face, pos, slots_d.pt[i], slots_d.q[i],
					slots_d.weight[i], slots_d.prompt_by[i], v, slots_d.mass,
					DeviceAdd {})};
				const double z_wall {Zones::face(zones, zone, Zones::zmin)
					== zone ? Zones::wrap_z(zones, zone, pos[2]) : pos[2]};
				if (self_sputter(e, z_wall)) state = 0;
			}

			// Absorb particles that ended up in a cell inside the wall, and
			// put them on the element bordering it whose plasma cell they
			// came from (the slot's zone and cell from before this step)
			int idx[3];
			if (!state && p.chi != nullptr
				&& in_wall_cell(zones, p.chi, zone, pos, idx))
			{
				state = 1;
				++tally.n[fate_chi];
				if (deposit)
				{
					const long e {Sputtering::add_into_wall_cell(dep, zone,
						idx[0], idx[1], idx[2], slots_d.zone[i],
						slots_d.xidx[i], slots_d.yidx[i], slots_d.pt[i],
						slots_d.q[i], slots_d.weight[i], slots_d.prompt_by[i],
						v, slots_d.mass, DeviceAdd {})};
					if (self_sputter(e, pos[2])) state = 0;
				}
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
	* @param dep_d Deposition tallies on this GPU
	*        (Sputtering::Deposition::to_device), or nullptr to not tally.
	* @param self_d Self-sputtering on this GPU
	*        (Sputtering::SelfSource::to_device), or nullptr for none.
	*        Needs dep_d.
	*/
	void check_bounds_gpu(Slots::SlotsDevice& slots_d,
		const Background::BackgroundDevice& bkg_d, const BoundaryParams& p,
		pcg32* rngs_d, unsigned long long* fates_d,
		const Sputtering::DepositionDevice* dep_d,
		const Sputtering::SelfDevice* self_d)
	{
		const int blockSize {256};
		const int gridSize {(slots_d.N + blockSize - 1) / blockSize};
		const bool deposit {dep_d != nullptr && dep_d->device_id >= 0};
		const bool self_on {deposit && self_d != nullptr
			&& self_d->device_id >= 0};
		check_bounds_kernel<<<gridSize, blockSize>>>(slots_d, bkg_d, p,
			rngs_d, fates_d, deposit ? dep_d->view
			: Sputtering::DepositionView {}, deposit, self_on ? self_d->view
			: Sputtering::SelfView {}, self_on, self_on ? self_d->avalanche
			: nullptr);

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
		fates.self_sputtered += static_cast<long long>(
			n[fate_self_sputtered]);
	}

}  // namespace Boundary
