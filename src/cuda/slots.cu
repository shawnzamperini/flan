/**
* @file  slots.cu
*
* @brief CUDA routines for filling particle slots with new particles and
*        checking if every slot is dead.
*
* GPU version of the slot routines in slots.cpp. Each time step, the
* indices of the dead slots are listed in increasing order
* (cub::DeviceSelect), and then one thread per dead slot creates a new
* particle there, as long as there are particles left to launch. New
* particles start according to the imp_*start_opt and imp_temp_start_opt
* options, or from the sputtered source (imp_source_opt = "sputter").
*/
#include <cstdio>

#include <cub/device/device_select.cuh>
#include <thrust/iterator/counting_iterator.h>

#include "background_device.h"
#include "errors.h"
#include "options_device.h"
#include "pcg32.h"
#include "slots.cuh"
#include "slots_device.h"
#include "sputtering_hd.h"
#include "zones_hd.h"

#include "indices.cuh"
#include "interpolate.cuh"


namespace Slots
{

	/**
	* Simple struct to hold the three components of a vector.
	*/
	struct Vec3d
	{
		double x;
		double y;
		double z;
	};

	/**
	* Picks the starting value of one coordinate (t, x, y or z) of a new
	* particle.
	*
	* start_opt_int is the control integer of the imp_*start_opt option
	* deciding where the particle starts:
	*   0: At start_val.
	*   1: Uniformly random between range_min and range_max.
	*   2: Uniformly random over the full simulation volume, between bkg_min
	*      and bkg_max.
	* Any other value returns 0.
	*
	* @param start_opt_int Control integer of the imp_*start_opt option.
	* @param start_val Starting value used when start_opt_int is 0.
	* @param range_min Lower end of the range used when start_opt_int is 1.
	* @param range_max Upper end of the range used when start_opt_int is 1.
	* @param bkg_min Minimum of the simulation volume, used when start_opt_int
	*        is 2.
	* @param bkg_max Maximum of the simulation volume, used when start_opt_int
	*        is 2.
	* @param rng Random number generator for this thread.
	* @return Starting value of the coordinate.
	*/
	__device__
	double get_birth_val_cuda(const int start_opt_int, const double start_val, 
		const double range_min, const double range_max, const double bkg_min, 
		const double bkg_max, pcg32& rng)
	{
		// Start at specific point
		double return_start_val {};
		if (start_opt_int == 0)
		{
			 return_start_val = start_val;
		}

		// Start between a user-specified range 
		else if (start_opt_int == 1)
		{
			return_start_val = range_min + (range_max - range_min) 
				* rng.next_double();
		}

		// Start between the full range of the simulation volume 
		else if (start_opt_int == 2)
		{
			return_start_val = bkg_min + (bkg_max - bkg_min) 
				* rng.next_double();
		}
		
		return return_start_val;
	}


	/**
	* Picks the starting Cartesian velocity of a new particle.
	*
	* Each component is sampled from a Maxwellian with zero mean
	* and sigma = sqrt(kT/m), giving an isotropic velocity distribution. The
	* temperature is imp_temp_start_val if imp_temp_start_opt is
	* "single_value", or the main ion temperature interpolated (trilinear
	* in space, linear in time) at the particle location if it is
	* "main_ion". The interpolated ion temperature is held at or above
	* min_ti, since it can go negative near the grid edges. The directed
	* starting velocity (imp_vX_start_val, imp_vY_start_val,
	* imp_vZ_start_val) is added to the sampled one.
	*
	* @param bkg_d Background plasma on the device, used for the ion
	*        temperature.
	* @param t Starting time of the particle (s).
	* @param x Starting computational x coordinate.
	* @param y Starting computational y coordinate.
	* @param z Starting computational z coordinate.
	* @param zone Zone the particle starts in.
	* @param tidx Index of the nearest background time frame.
	* @param xidx x index of the grid cell containing the particle.
	* @param yidx y index of the grid cell containing the particle.
	* @param zidx z index of the grid cell containing the particle.
	* @param rng Random number generator for this thread.
	* @param opts_d Simulation options on the device.
	* @return Starting velocity {vX, vY, vZ} in Cartesian space (m/s).
	*/
	__device__
	Vec3d get_birth_vXYZ(const Background::BackgroundDevice bkg_d, 
		const double t, const double x, const double y, const double z, 
		const int zone, const int tidx, const int xidx, const int yidx,
		const int zidx, pcg32& rng, const Options::OptionsDevice* opts_d)
	{
		// Start at input temperature value
		double start_temp {};
		if (opts_d->imp_temp_start_opt_int == 0)
		{
			start_temp = opts_d->imp_temp_start_val;
		}

		// Start at main ion temperature, interpolating for it
		// Get interpolation stencil at this location for 4D arrays
		else if (opts_d->imp_temp_start_opt_int == 1)
		{
			Interpolate::InterpolationStencil stencil_4d
				{Interpolate::build_stencil_4d(bkg_d, zone, tidx, xidx, yidx, zidx, 
				t, x, y, z)};

			start_temp = Interpolate::interp_plasma(bkg_d, bkg_d.ti,
				stencil_4d, zone, t, x, y, z);

			// Limit to minimum value
			start_temp = fmax(start_temp, opts_d->min_ti);
		}

		// Sample from a Maxwellian with sigma = sqrt(kT/m) and mean = 0 for an
		// isotropic velocity distribution. T [eV], m [kg]
		const double sigma {sqrt(start_temp * Constants::ev_to_j 
			/ (opts_d->imp_mass_amu * Constants::amu_to_kg))};  // m/s

		double vX {rng.normal(0.0, sigma)};
		double vY {rng.normal(0.0, sigma)};
		double vZ {rng.normal(0.0, sigma)};

		// Plus the directed starting velocity (imp_v*_start_val, 0 by
		// default), e.g., for a beam or a case where every particle must
		// start the same way
		vX += opts_d->imp_vX_start_val;
		vY += opts_d->imp_vY_start_val;
		vZ += opts_d->imp_vZ_start_val;

		return {vX, vY, vZ};
	}

	
	/**
	* Creates a new particle based on the input options.
	*
	* The particle starts in zone imp_zone_start_val (as
	* Slots::make_new_particle on the CPU), and its starting t, x, y and z
	* come from get_birth_val_cuda in that zone's coordinates, the indices
	* of the nearest time frame and the containing grid cell are found, and
	* the starting Cartesian velocity comes from get_birth_vXYZ. The
	* computational velocities are left at zero, since they are set on the
	* first Boris update. The weight is imp_source_scale_fact / imp_num,
	* the atoms/s the particle stands for, and the charge is init_charge.
	*
	* @param rng_d Random number generator of the slot the particle goes
	*        into.
	* @param opts_d Simulation options on the device.
	* @param bkg_d Background plasma on the device.
	* @return The new particle.
	*/
	__device__
	ParticleInitDevice make_new_particle_cuda(pcg32& rng, 
		const Options::OptionsDevice* opts_d, 
		const Background::BackgroundDevice& bkg_d)
	{

		ParticleInitDevice p;

		// The zone the particle is born in (imp_zone_start_val, numbered
		// from 1 in the input). The starting position is in that zone's
		// coordinates.
		const Zones::ZonesView& zones {bkg_d.zones};
		p.zone = opts_d->zone_start;
		const double* gx {Zones::zone_edges(zones, p.zone, 0)};
		const double* gy {Zones::zone_edges(zones, p.zone, 1)};
		const double* gz {Zones::zone_edges(zones, p.zone, 2)};

		// Create a new particle based on input options. The full range is
		// the zone's outer cell edges.
		p.t = get_birth_val_cuda(opts_d->tstart_opt_int, opts_d->tstart_val, 
			opts_d->trange_min, opts_d->trange_max, bkg_d.t_min,
			bkg_d.t_max, rng);
		p.x = get_birth_val_cuda(opts_d->xstart_opt_int, opts_d->xstart_val, 
			opts_d->xrange_min, opts_d->xrange_max, gx[0],
			gx[zones.n[0][p.zone]], rng);
		p.y = get_birth_val_cuda(opts_d->ystart_opt_int, opts_d->ystart_val, 
			opts_d->yrange_min, opts_d->yrange_max, gy[0],
			gy[zones.n[1][p.zone]], rng);
		p.z = get_birth_val_cuda(opts_d->zstart_opt_int, opts_d->zstart_val, 
			opts_d->zrange_min, opts_d->zrange_max, gz[0],
			gz[zones.n[2][p.zone]], rng);

		// Get indices, same process as in impurity_transport.cu
		p.tidx = Indices::get_nearest_index_cuda(bkg_d.times, bkg_d.tdim,
			p.t);
		p.xidx = Zones::cell_index(zones, p.zone, 0, p.x);
		p.yidx = Zones::cell_index(zones, p.zone, 1, p.y);
		p.zidx = Zones::cell_index(zones, p.zone, 2, p.z);

		// Get starting vX, vY, vZ into struct
		Vec3d vXYZ {get_birth_vXYZ(bkg_d, p.t, p.x, p.y, p.z, p.zone, p.tidx,
			p.xidx, p.yidx, p.zidx, rng, opts_d)};
		p.vX = vXYZ.x;
		p.vY = vXYZ.y;
		p.vZ = vXYZ.z;

		// These get set on the first call to Boris::update_velocity 
		p.vx = 0.0;
		p.vy = 0.0;
		p.vz = 0.0;

		// The weight is the atoms/s this particle stands for: the source,
		// imp_source_scale_fact, shared between the imp_num particles
		p.weight = opts_d->imp_source_scale_fact / opts_d->imp_num;
		p.q = opts_d->init_charge;

		return p;
	}

	
	/**
	* Starts a sputtered atom from the source, as make_new_particle does on
	* the CPU (slots.cpp): where, when and how fast from
	* Sputtering::sample_birth (the same routine, with the same draws),
	* neutral, and standing for source_rate / imp_num atoms/s.
	*
	* @param rng Random number generator of the slot.
	* @param source_d The source on the device.
	* @param bkg_d Background plasma on the device.
	* @return The new particle.
	*/
	__device__
	ParticleInitDevice make_sputtered_particle_cuda(pcg32& rng,
		const Sputtering::SourceDevice& source_d,
		const Background::BackgroundDevice& bkg_d)
	{
		const Sputtering::Birth b {Sputtering::sample_birth(source_d.view,
			rng)};
		ParticleInitDevice p;
		p.zone = b.zone;
		p.t = b.t;
		p.x = b.x;
		p.y = b.y;
		p.z = b.z;
		p.tidx = Indices::get_nearest_index_cuda(bkg_d.times, bkg_d.tdim,
			p.t);
		p.xidx = Zones::cell_index(bkg_d.zones, p.zone, 0, p.x);
		p.yidx = Zones::cell_index(bkg_d.zones, p.zone, 1, p.y);
		p.zidx = Zones::cell_index(bkg_d.zones, p.zone, 2, p.z);
		p.vX = b.vX;
		p.vY = b.vY;
		p.vZ = b.vZ;
		p.vx = 0.0;
		p.vy = 0.0;
		p.vz = 0.0;
		p.weight = source_d.weight;
		p.q = 0;
		return p;
	}

	/**
	* Predicate for cub::DeviceSelect::If: true for the index of a dead
	* slot (state not 0).
	*/
	struct IsDead
	{
		const int* state;	///< State of each slot, 0 if alive

		/**
		* @param i Slot index.
		* @return true if slot i is dead.
		*/
		__device__ bool operator()(const int i) const {return state[i] != 0;}
	};

	/// Selection predicate for cub::DeviceSelect::If: whether a slot holds
	/// a live particle (state 0), for the list of live slots
	struct IsAlive
	{
		const int* state;	///< State of each slot, 0 if alive

		/**
		* @param i Slot index.
		* @return true if slot i is alive.
		*/
		__device__ bool operator()(const int i) const {return state[i] == 0;}
	};


	/**
	* Kernel that fills dead slots with new particles.
	*
	* One thread handles one entry of the dead slot list from
	* fill_slots_gpu, so only the first min(num_dead, rem_parts) dead
	* slots are filled. Each creates a new particle with
	* make_new_particle_cuda and writes it into the slot, marking it alive.
	* The cell indices are set to 0 here and assigned afterward by
	* find_containing_cell_gpu.
	*
	* slots_d is passed by value each time since the copy time is much less
	* than the kernel launch time, and it makes for clearer code.
	*
	* @param slots_d Particle slots on the device. The dead slots are
	*        overwritten in place.
	* @param bkg_d Background plasma on the device.
	* @param rem_parts Number of particles left to launch.
	* @param rngs_d Random number generators on the device.
	* @param opts_d Simulation options on the device.
	* @param source Sputtered source on the device, or a copy with
	*        device_id -1 to start particles from the imp_*start_opt
	*        options.
	*/
	__global__ void fill_slots_kernel(SlotsDevice slots_d, 
		Background::BackgroundDevice bkg_d, int rem_parts, pcg32* rngs_d, 
		const Options::OptionsDevice* opts_d,
		const Sputtering::SourceDevice source)
	{
		// Global index
		int idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx >= *slots_d.num_dead) return;
		if (idx >= rem_parts) return;

		// Pick slot index only from those that have been found to be dead.
		int i = slots_d.dead_indices[idx];

		// Create a new particle (device-side initializer)
		ParticleInitDevice p = source.device_id >= 0
			? make_sputtered_particle_cuda(rngs_d[i], source, bkg_d)
			: make_new_particle_cuda(rngs_d[i], opts_d, bkg_d);

		// Write particle data
		slots_d.pt[i] = 0.0;
		slots_d.t[i]  = p.t;

		slots_d.x[i]  = p.x;
		slots_d.y[i]  = p.y;
		slots_d.z[i]  = p.z;

		slots_d.vx[i] = p.vx;
		slots_d.vy[i] = p.vy;
		slots_d.vz[i] = p.vz;

		slots_d.vX[i] = p.vX;
		slots_d.vY[i] = p.vY;
		slots_d.vZ[i] = p.vZ;
		
		slots_d.weight[i] = p.weight;

		slots_d.q[i] = p.q;
		slots_d.zone[i] = p.zone;

		// The time frame and cell of the birth, as fill_slots_cpu (the
		// main loop finds them again with find_containing_cell_gpu)
		slots_d.tidx[i] = p.tidx;
		slots_d.xidx[i] = p.xidx;
		slots_d.yidx[i] = p.yidx;
		slots_d.zidx[i] = p.zidx;

		// Not yet ionized (no prompt window) and not self-sputtered, as
		// fill_slots_cpu
		slots_d.prompt_by[i] = -1.0;
		slots_d.generation[i] = 0;

		// Mark slot alive
		slots_d.state[i] = 0;
	}

	/**
	* Kernel that sets all_dead to false if any slot holds an alive
	* particle.
	*
	* all_dead must be set to true before the launch. One thread checks
	* one slot.
	*
	* @param slots_d Particle slots on the device.
	*/
	__global__ void all_dead_kernel(SlotsDevice slots_d)
	{
		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= slots_d.N) return;

		// If a single thread has an alive particle (state = 0), then all_dead
		// must be false. This is technically a race condition, but we allow
		// it since every thread that writes is writing the same value (false).
		if (slots_d.state[i] == 0) *slots_d.all_dead = false;
	}

	/**
	* Replaces dead particles with new ones, as long as there are particles
	* left to launch.
	*
	* Lists the dead slot indices in increasing order with
	* cub::DeviceSelect::If, copies the number of dead slots back to the
	* host (which waits for the selection to finish), and then launches
	* fill_slots_kernel with one thread per dead slot, in blocks of 256
	* threads.
	*
	* The order matters when there are more dead slots than particles left
	* to launch: only the first rem_parts get a particle, and each slot has
	* its own random number generator. An order that is the same every run
	* (unlike one from atomics) makes a GPU run reproducible, so changes to
	* the code can be checked bit for bit. The selection's work space is
	* allocated here the first time and kept in slots_d. The fill kernel launch is
	* asynchronous, except in DEBUG builds, where it waits for the kernel to
	* finish and prints any error. The number of filled slots is
	* min(rem_parts, num_dead), so the counts are updated on the host
	* without waiting for it.
	*
	* @param slots_d Particle slots on the device.
	* @param rem_parts Number of particles left to launch. Decreased by the
	*        number of slots filled.
	* @param[out] alive_slots Number of alive slots after filling.
	* @param rngs_d Random number generators on the device.
	* @param opts_d Simulation options on the device.
	* @param bkg_d Background plasma on the device.
	* @param source_d Sputtered source on the device to start particles
	*        from, or nullptr to start them from the imp_*start_opt options.
	*/
	void fill_slots_gpu(SlotsDevice& slots_d, int& rem_parts, int& alive_slots,
		pcg32* rngs_d, Options::OptionsDevice* opts_d,
		const Background::BackgroundDevice& bkg_d,
		const Sputtering::SourceDevice* source_d)
	{

		// Each slot is assigned to a specific GPU where its data lies
		cudaSetDevice(slots_d.device_id);

		// List the indices of the dead slots, in increasing order, in
		// dead_indices, and their number in num_dead. The first call (with
		// no work space) only gives the size of the work space needed.
		const thrust::counting_iterator<int> slot_index(0);
		const IsDead is_dead {slots_d.state};
		std::size_t temp_bytes {};
		cub::DeviceSelect::If(nullptr, temp_bytes, slot_index,
			slots_d.dead_indices, slots_d.num_dead, slots_d.N, is_dead);
		if (temp_bytes > slots_d.select_temp_bytes)
		{
			cudaFree(slots_d.select_temp);
			Errors::check_cuda(cudaMalloc(&slots_d.select_temp, temp_bytes),
				"fill_slots_gpu: allocating the work space");
			slots_d.select_temp_bytes = temp_bytes;
		}
		cub::DeviceSelect::If(slots_d.select_temp, slots_d.select_temp_bytes,
			slot_index, slots_d.dead_indices, slots_d.num_dead, slots_d.N,
			is_dead);

		// Read back num_dead
		int num_dead = 0;
		cudaMemcpy(&num_dead, slots_d.num_dead, sizeof(*slots_d.num_dead), 
			cudaMemcpyDeviceToHost);

		int rem_initial {rem_parts};

		// Call GPU kernel to fill slots.
		int fillBlock = 256;
		int fillGrid  = (num_dead + fillBlock - 1) / fillBlock;
		const Sputtering::SourceDevice source {source_d != nullptr
			? *source_d : Sputtering::SourceDevice {}};
		fill_slots_kernel<<<fillGrid, fillBlock>>>(slots_d, bkg_d, rem_parts,
			rngs_d, opts_d, source);

#ifdef DEBUG
		// Check for errors. This also catches errors from compact_dead_slots.
		cudaError_t err {cudaDeviceSynchronize()};
		if (err != cudaSuccess)
			printf("fill_slots_kernel error: %s\n", cudaGetErrorString(err));
#endif

		// Retrieve how many particles were actually filled
		int filled {std::min(rem_initial, num_dead)};

		// Retrieve number of alive slots
		alive_slots = slots_d.N - num_dead + filled;

		// Update remaining particles
		rem_parts = rem_initial - filled;

		// List the live slots, in increasing order, for the kernels until
		// the next fill (Slots::slot_of_thread). Only check_bounds_kernel
		// and step_kernel kill particles, and they run before the next
		// fill, so alive_slots is the length of the list and the launches
		// are sized from it without waiting for the selection.
		if (slots_d.live_indices == nullptr)
		{
			Errors::check_cuda(cudaMalloc(&slots_d.live_indices,
				slots_d.N * sizeof(int)),
				"fill_slots_gpu: allocating the live list");
			Errors::check_cuda(cudaMalloc(&slots_d.num_live_d, sizeof(int)),
				"fill_slots_gpu: allocating the live count");
		}
		const IsAlive is_alive {slots_d.state};
		std::size_t live_temp_bytes {};
		cub::DeviceSelect::If(nullptr, live_temp_bytes, slot_index,
			slots_d.live_indices, slots_d.num_live_d, slots_d.N, is_alive);
		if (live_temp_bytes > slots_d.select_temp_bytes)
		{
			cudaFree(slots_d.select_temp);
			Errors::check_cuda(cudaMalloc(&slots_d.select_temp,
				live_temp_bytes), "fill_slots_gpu: allocating the work space");
			slots_d.select_temp_bytes = live_temp_bytes;
		}
		cub::DeviceSelect::If(slots_d.select_temp, slots_d.select_temp_bytes,
			slot_index, slots_d.live_indices, slots_d.num_live_d, slots_d.N,
			is_alive);
		slots_d.num_live = alive_slots;

#ifdef DEBUG
		// The host's count must be the selection's
		int num_live_check {};
		cudaMemcpy(&num_live_check, slots_d.num_live_d, sizeof(int),
			cudaMemcpyDeviceToHost);
		if (num_live_check != alive_slots)
		{
			printf("fill_slots_gpu: %d live slots listed, %d expected\n",
				num_live_check, alive_slots);
		}
#endif
	}


	/**
	* Checks if every slot holds a dead particle.
	*
	* Sets all_dead to true on the device, launches all_dead_kernel (blocks
	* of 256 threads, one thread per slot) and copies the result back,
	* which waits for the kernel to finish. In DEBUG builds, any kernel
	* error is printed.
	*
	* @param slots_d Particle slots on the device.
	* @return true if every slot is dead, false if not.
	*/
	bool all_dead_gpu(SlotsDevice& slots_d)
	{

		// Each slot is assigned to a specific GPU where its data lies
		cudaSetDevice(slots_d.device_id);

		// Must set all_dead to true each time. Think of it like a question
		// that the device responds to.
		//   Host: Everything is dead? (all_dead = true)
		//   Device: No (all_dead = false)
		//   Device: Yes (all_dead untouched, stays true)
		// As we come back again during the next time step, we have to ask the
		// question again.
		bool init = true;
		cudaMemcpy(slots_d.all_dead, &init, sizeof(bool), 
			cudaMemcpyHostToDevice);

		// Call GPU kernel
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;
		all_dead_kernel<<<gridSize, blockSize>>>(slots_d);

#ifdef DEBUG
		// Check for errors
		cudaError_t err {cudaDeviceSynchronize()};
		if (err != cudaSuccess)
			printf("all_dead_kernel error: %s\n", cudaGetErrorString(err));
#endif

		// Copy all_dead back to host and return
		bool all_dead;
		cudaMemcpy(&all_dead, slots_d.all_dead, sizeof(bool), 
			cudaMemcpyDeviceToHost);

		return all_dead;
	}
}
