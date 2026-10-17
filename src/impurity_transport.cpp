/**
* @file impurity_transport.cpp
* @brief Contains main impurity following loop and wrappers for calling the
*        various steps in following particles.
*/

#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <omp.h>
#include <thread>
#include <tuple>
#include <vector>

#include "background.h"
#include "boris.h"
#include "boundary.h"
#include "collisions.h"
#include "constants.h"
#include "flan_types.h"
#include "impurity_stats.h"
#include "impurity_transport.h"
#include "indices.h"
#include "openadas.h"
#include "options.h"
#include "options_device.h"
#include "particle_progress.h"
#include "pcg32.h"
#include "slots.h"
#include "timer.h"
#include "utilities.h"
#include "vectors.h"

#ifdef USE_CUDA
#include <cuda_runtime.h>
#include "boris.cuh"
#include "boundary.cuh"
#include "collisions.cuh"
#include "init_rngs.cuh"
#include "impurity_stats.cuh"
#include "impurity_transport.cuh"
#include "openadas.cuh"
#include "slots.cuh"
#endif


namespace ImpurityTransport
{
	/**
	* Finds the background time index and spatial cell indices for every
	* particle in slots (CPU implementation).
	*
	* Time and space are handled differently. The background times are
	* stored at "time centers" rather than on a grid, so the time index is
	* the index of the background time nearest to the particle's time. The
	* spatial coordinates are defined by grid edges, so each spatial index
	* is the index of the cell that contains the particle's x, y and z
	* coordinates. Particles outside the grid are assigned the first or
	* last cell in that dimension.
	*
	* The indices are written back into slots with set_tidx, set_xidx,
	* set_yidx and set_zidx. The loop over particles is parallelized with
	* OpenMP. Every slot is updated, including slots holding dead
	* particles.
	*
	* @param slots Particle slots whose positions are read and whose time
	*              and cell indices are updated.
	* @param bkg Background plasma that provides the time array and the
	*            x, y and z grid edges.
	*/
	void find_containing_cell_cpu(Slots::Slots& slots, 
		const Background::Background& bkg)
	{
		// Variables for each grid to avoid repeatedly calling in the loop
		const auto& times = bkg.get_times();
		const auto& grid_x = bkg.get_grid_x();
		const auto& grid_y = bkg.get_grid_y();
		const auto& grid_z = bkg.get_grid_z();

		// Loop through each slot
		#pragma omp parallel for
		for (int i = 0; i < slots.N(); i++)
		{
			// We don't have a "time grid" since time is stored at 
			// "time centers", if that helps think about it. So we just find
			// the nearest index for it. The spatial coordinates have a grid
			// so we use that.
			int tidx {Indices::get_nearest_index_cpu(times, slots.t()[i])};

			// Find nearest cell using the grid (unlike time, which uses
			// cell "centers"). 
			int xidx {Indices::get_nearest_cell_index_cpu(grid_x, slots.x()[i])};
			int yidx {Indices::get_nearest_cell_index_cpu(grid_y, slots.y()[i])};
			int zidx {Indices::get_nearest_cell_index_cpu(grid_z, slots.z()[i])};

			slots.set_tidx(i, tidx);
			slots.set_xidx(i, xidx);
			slots.set_yidx(i, yidx);
			slots.set_zidx(i, zidx);
		}
	}

	
	/**
	* Advances every alive particle in slots by one time step (CPU
	* implementation).
	*
	* Each alive particle's time since launch (pt) and background time (t)
	* are increased by dt, and its x, y and z positions are advanced using
	* its current velocity. If the particle's time since launch exceeds the
	* imp_max_time option, the particle is marked dead and max_time_ctr is
	* incremented. The iteration counter of each stepped particle is also
	* incremented. Dead particles are skipped. The loop over particles is
	* parallelized with OpenMP.
	*
	* @param slots Particle slots to advance.
	* @param bkg Background plasma (currently unused).
	* @param opts Simulation options, used for imp_max_time.
	* @param dt Time step (s).
	* @param max_time_ctr Counter of particles that exceeded the maximum
	*                     allowed time. Incremented atomically.
	*/
	void step_cpu(Slots::Slots& slots, const Background::Background& bkg, 
		const Options::Options& opts, const double dt, int& max_time_ctr)
	{

		#pragma omp parallel for
		for (int i=0; i < slots.N(); ++i)
		{
			// Don't step dead particles
			if (slots.state()[i]) continue;

			/*
			if (i == 0)
			{
				#pragma omp critical
				{
					std::cout << "======================\n";
					std::cout << "i = " << i << '\n';
					std::cout << "q = " << slots.q()[i] << '\n';
					std::cout << "tidx = " << slots.tidx()[i] << '\n';
					std::cout << "xidx = " << slots.xidx()[i] << '\n';
					std::cout << "yidx = " << slots.yidx()[i] << '\n';
					std::cout << "zidx = " << slots.zidx()[i] << '\n';
					std::cout << "t = " << slots.t()[i] << '\n';
					std::cout << "x = " << slots.x()[i] << '\n';
					std::cout << "y = " << slots.y()[i] << '\n';
					std::cout << "z = " << slots.z()[i] << '\n';
					std::cout << "vx = " << slots.vx()[i] << '\n';
					std::cout << "vy = " << slots.vy()[i] << '\n';
					std::cout << "vz = " << slots.vz()[i] << '\n';
					std::cout << "vX = " << slots.vX()[i] << '\n';
					std::cout << "vY = " << slots.vY()[i] << '\n';
					std::cout << "vZ = " << slots.vZ()[i] << '\n';
				}
			}
			*/

			// Update particle time and time in the background
			slots.set_pt(i, slots.pt()[i] + dt);
			slots.set_t(i, slots.t()[i] + dt);

			// Update curvilinear position
			slots.set_x(i, slots.x()[i] + slots.vx()[i] * dt);
			slots.set_y(i, slots.y()[i] + slots.vy()[i] * dt);
			slots.set_z(i, slots.z()[i] + slots.vz()[i] * dt);

			// Check if maximum allowed time has been exceeded, setting to
			// dead if so and incrementing counter
			if (slots.pt()[i] > opts.imp_max_time())
			{
				slots.set_state(i, 1);

				#pragma omp atomic
				max_time_ctr++;
			}

			// Update number of iterations (mainly for debugging)
			slots.set_iters(i, slots.iters()[i] + 1);
		}
	}

	
	/**
	* Updates particle charge states from ionization and recombination
	* probabilities (CPU implementation).
	*
	* This is a thin wrapper around OpenADAS::ioniz_recomb, which uses the
	* OpenADAS rate coefficients and the local background plasma to decide
	* whether each particle ionizes or recombines during the time step. The
	* particle charges in slots are modified in place.
	*
	* @param slots Particle slots whose charge states are updated.
	* @param bkg Background plasma that provides the local density and
	*            temperature.
	* @param oa_ioniz OpenADAS ionization rate coefficients (SCD).
	* @param oa_recomb OpenADAS recombination rate coefficients (ACD).
	* @param opts Simulation options (currently unused).
	* @param dt Time step (s).
	* @param ioniz_warnings Counter of particles whose ionization
	*                       probability exceeded 1.
	* @param recomb_warnings Counter of particles whose recombination
	*                        probability exceeded 1.
	* @param rngs Random number generators, one per OpenMP thread.
	*/
	void ioniz_recomb_cpu(Slots::Slots& slots, 
		const Background::Background& bkg, const OpenADAS::OpenADAS& oa_ioniz,
		const OpenADAS::OpenADAS& oa_recomb, const Options::Options& opts, 
		const double dt, int& ioniz_warnings, int& recomb_warnings, 
		std::vector<pcg32>& rngs)
	{

		// Updates particle charge within
		OpenADAS::ioniz_recomb(slots, bkg, oa_ioniz, oa_recomb, dt, 
			ioniz_warnings, recomb_warnings, rngs);
	}


	/**
	* Applies Coulomb collisions with the background ions and electrons to
	* every particle in slots (CPU implementation).
	*
	* Particle velocities are updated in place with the Nanbu collision
	* model. Collisions with the background ions are always applied.
	* Collisions with the electrons are applied afterwards, except for the
	* friction_force test case, which only considers ion collisions so the
	* result can be compared against the expected flow.
	*
	* @param slots Particle slots whose velocities are updated.
	* @param bkg Background plasma that provides the ion and electron
	*            properties.
	* @param dt Time step (s).
	* @param opts Simulation options, used to check the test_opt option.
	* @param imp_stats Impurity statistics passed through to the collision
	*                  model.
	* @param rngs Random number generators, one per OpenMP thread.
	*/
	void collision_cpu(Slots::Slots& slots, const Background::Background& bkg, 
		const double dt, const Options::Options& opts, 
		Impurity::Statistics& imp_stats, std::vector<pcg32>& rngs)
	{
		// Update impurity velocity based on Nanbu collision model. Impurity
		// is modified within function. First call is for ions (the false) and
		// second call is for electrons (the true).
		Collisions::nanbu_coll(slots, bkg, opts, false, dt, imp_stats, rngs);

		// friction_force test case only considers ion collisions to compare
		// against expected flow
		if (opts.test_opt_int() != 5)
		{
			Collisions::nanbu_coll(slots, bkg, opts, true, dt, imp_stats, rngs);
		}
	}


	// ------------------
	// Wrapper functions
	// ------------------

	/**
	* Finds the background time index and spatial cell indices for every
	* particle, choosing the CPU or GPU implementation.
	*
	* If GPU acceleration is on (use_gpu), find_containing_cell_gpu is called
	* on the device-side slots. Otherwise find_containing_cell_cpu is called
	* on the host-side slots. The GPU branch only exists in CUDA builds.
	*
	* @param slots Host-side particle slots (used on the CPU path).
	* @param slots_d Device-side particle slots (used on the GPU path).
	* @param bkg Host-side background plasma (used on the CPU path).
	* @param bkg_d Device-side background plasma (used on the GPU path).
	* @param opts Simulation options, used to check whether the GPU is on.
	*/
	void find_containing_cell_wrapper(Slots::Slots& slots, 
		Slots::SlotsDevice& slots_d, const Background::Background& bkg, 
		const Background::BackgroundDevice& bkg_d, const Options::Options& opts)
	{

#ifdef USE_CUDA
		if (opts.use_gpu_int() > 0)
		{
			// Defined in impurity_transport.cu
			find_containing_cell_gpu(slots_d, bkg_d);
			return;
		}
#endif

		find_containing_cell_cpu(slots, bkg);
	}


	/**
	* Records the current particle state in the impurity statistics arrays,
	* choosing the CPU or GPU implementation.
	*
	* If GPU acceleration is on (use_gpu), ImpurityStats::record_stats_gpu
	* records into the device-side statistics. Otherwise record_stats_cpu
	* records into the host-side statistics. The GPU branch only exists in
	* CUDA builds.
	*
	* @param imp_stats Host-side impurity statistics (used on the CPU path).
	* @param imp_stats_d Device-side impurity statistics (used on the GPU
	*                    path).
	* @param slots Host-side particle slots (used on the CPU path).
	* @param slots_d Device-side particle slots (used on the GPU path).
	* @param opts Simulation options, used to check whether the GPU is on.
	* @param imp_time_step Impurity time step (s).
	*/
	void record_stats_wrapper(Impurity::Statistics& imp_stats, 
		ImpurityStats::StatisticsDevice& imp_stats_d, 
		const Slots::Slots& slots, const Slots::SlotsDevice& slots_d, 
		const Options::Options& opts, const double imp_time_step)
	{

#ifdef USE_CUDA
		if (opts.use_gpu_int() > 0)
		{
			// Defined in cuda/impurity_stats.cu
			ImpurityStats::record_stats_gpu(imp_stats_d, slots_d, imp_time_step);
			return;
		}
#endif

		// Defined in impurity_stats.cpp
		record_stats_cpu(imp_stats, slots, opts, imp_time_step);
	}


	/**
	* Replaces dead particles in slots with newly launched particles,
	* choosing the CPU or GPU implementation.
	*
	* Dead slots are revived with new particles until either every slot is
	* occupied or there are no particles left to launch. The number of
	* particles launched is subtracted from rem_parts. If GPU acceleration
	* is on (use_gpu), Slots::fill_slots_gpu is called. Otherwise
	* Slots::fill_slots_cpu is called. The GPU branch only exists in CUDA
	* builds.
	*
	* @param slots Host-side particle slots (used on the CPU path).
	* @param slots_d Device-side particle slots (used on the GPU path).
	* @param rem_parts Number of particles still waiting to be launched.
	*                  Reduced by the number of particles launched.
	* @param bkg Host-side background plasma (used on the CPU path).
	* @param opts Simulation options.
	* @param alive_slots Set to the number of alive particles in slots,
	*                    used for progress printouts.
	* @param rngs Host-side random number generators, one per OpenMP
	*             thread (used on the CPU path).
	* @param rngs_d Device-side random number generators (used on the GPU
	*               path).
	* @param opts_d Device-side simulation options (used on the GPU path).
	* @param bkg_d Device-side background plasma (used on the GPU path).
	*/
	void fill_slots_wrapper(Slots::Slots& slots, Slots::SlotsDevice& slots_d,
		int& rem_parts, const Background::Background& bkg, 
		Options::Options& opts, int& alive_slots, std::vector<pcg32>& rngs,
		pcg32* rngs_d, Options::OptionsDevice* opts_d, 
		const Background::BackgroundDevice& bkg_d)
	{

#ifdef USE_CUDA
		if (opts.use_gpu_int() > 0) 
		{
			// Defined in cuda/slots.cu
			Slots::fill_slots_gpu(slots_d, rem_parts, alive_slots, rngs_d, 
				opts_d, bkg_d);
			return;
		}
#endif

		// Defined in slots.cpp
		Slots::fill_slots_cpu(slots, rem_parts, alive_slots, bkg, opts, rngs);
	}


	/**
	* Updates particle velocities with the Boris algorithm, choosing the CPU
	* or GPU implementation.
	*
	* If GPU acceleration is on (use_gpu), Boris::update_velocity_gpu is
	* called. Otherwise Boris::update_velocity_cpu is called. The GPU branch
	* only exists in CUDA builds. A negative dt pushes the velocities
	* backward in time, which main_loop uses to offset the starting
	* velocities by half a time step.
	*
	* @param slots Host-side particle slots (used on the CPU path).
	* @param slots_d Device-side particle slots (used on the GPU path).
	* @param bkg Host-side background plasma (used on the CPU path).
	* @param bkg_d Device-side background plasma (used on the GPU path).
	* @param opts Simulation options.
	* @param dt Time step to advance the velocities by (s).
	*/
	void boris_wrapper(Slots::Slots& slots, 
		Slots::SlotsDevice& slots_d, const Background::Background& bkg, 
		const Background::BackgroundDevice& bkg_d, const Options::Options& opts,
		const double dt)
	{

#ifdef USE_CUDA
		if (opts.use_gpu_int() > 0) 
		{
			// Defined in cuda/boris.cu
			Boris::update_velocity_gpu(slots_d, bkg_d, dt);
			return;
		}
#endif

		Boris::update_velocity_cpu(slots, bkg, opts, dt);
	}


	/**
	* Checks whether the simulation is finished, choosing the CPU or GPU
	* implementation.
	*
	* The simulation is finished when no particles are waiting to be
	* launched and every particle in slots is dead. If rem_parts is
	* positive, false is returned without checking slots. If GPU
	* acceleration is on (use_gpu), Slots::all_dead_gpu is called, which
	* copies data from the GPU and so implicitly synchronizes it. The GPU
	* branch only exists in CUDA builds.
	*
	* @param slots Host-side particle slots (used on the CPU path).
	* @param slots_d Device-side particle slots (used on the GPU path).
	* @param opts Simulation options, used to check whether the GPU is on.
	* @param rem_parts Number of particles still waiting to be launched.
	* @return True if no particles remain to be launched and every slot is
	*         dead, false otherwise.
	*/
	bool all_dead_wrapper(Slots::Slots& slots, Slots::SlotsDevice& slots_d,
		const Options::Options& opts, int rem_parts)
	{
		// Still particles remaining to be placed into slots, not done yet.
		if (rem_parts > 0) return false;

#ifdef USE_CUDA
		if (opts.use_gpu_int() > 0) 
		{
			// Need to call a kernel here since the memory lives on the GPU
			return Slots::all_dead_gpu(slots_d);
		}
#endif

		// For CPU this is just a one-liner. If every entry is not equal
		// to zero (i.e., alive) then they're all dead.
		return std::all_of(slots.state().begin(), slots.state().end(),
			[](int s){return s != 0;});

	}


	/**
	* Checks whether particles have left the time or spatial bounds of the
	* background and handles them, choosing the CPU or GPU implementation.
	*
	* What happens to a particle that crosses a boundary is set by the
	* boundary type options (tbound_type, min_xbound_type, max_xbound_type,
	* etc.). If GPU acceleration is on (use_gpu), Boundary::check_bounds_gpu
	* is called with the boundary options unpacked from opts. Otherwise
	* Boundary::check_bounds_cpu is called. The GPU branch only exists in
	* CUDA builds.
	*
	* @param slots Host-side particle slots (used on the CPU path).
	* @param slots_d Device-side particle slots (used on the GPU path).
	* @param bkg Host-side background plasma (used on the CPU path).
	* @param bkg_d Device-side background plasma (used on the GPU path).
	* @param opts Simulation options, which provide the boundary types and
	*             buffers.
	* @param rngs_d Device-side random number generators (used on the GPU
	*               path).
	*/
	void check_bounds_wrapper(Slots::Slots& slots, Slots::SlotsDevice& slots_d,
		const Background::Background& bkg, 
		const Background::BackgroundDevice& bkg_d, 
		const Options::Options& opts, std::vector<pcg32>& rngs, pcg32* rngs_d)
	{

#ifdef USE_CUDA
		if (opts.use_gpu_int() > 0) 
		{
			// Defined in cuda/boundary.cu
			Boundary::check_bounds_gpu(slots_d, bkg_d, 
				opts.tbound_type_int(), 
				opts.min_xbound_type_int(), opts.max_xbound_type_int(),
				opts.min_ybound_type_int(), opts.max_ybound_type_int(),
				opts.min_zbound_type_int(), opts.max_zbound_type_int(),
				opts.imp_xbound_buffer(), opts.imp_ybound_buffer(),
				opts.imp_zbound_buffer(), opts.lcfs_x(), rngs_d);
			return;
		}
#endif

		// Defined in boundary.cpp
		Boundary::check_bounds_cpu(slots, bkg, opts, rngs);
	}

	
	/**
	* Advances every alive particle by one impurity time step, choosing the
	* CPU or GPU implementation.
	*
	* The time step is the imp_time_step option. If GPU acceleration is on
	* (use_gpu), step_gpu is called. Otherwise step_cpu is called. The GPU
	* branch only exists in CUDA builds. Both paths kill particles that
	* have been followed for longer than the imp_max_time option and count
	* them. The CPU path adds to max_time_ctr directly. The GPU path counts
	* them in warnings_d.max_time on the device, which is added to
	* max_time_ctr with OpenADAS::add_warnings_to_host after the run.
	*
	* @param slots Host-side particle slots (used on the CPU path).
	* @param slots_d Device-side particle slots (used on the GPU path).
	* @param bkg Host-side background plasma (used on the CPU path).
	* @param bkg_d Device-side background plasma (currently unused).
	* @param opts Simulation options, which provide the time step.
	* @param opts_d Device-side simulation options (used on the GPU path).
	* @param warnings_d Device-side warning counters for this GPU (used on
	*                   the GPU path).
	* @param max_time_ctr Counter of particles that exceeded the maximum
	*                     allowed time (CPU path only).
	*/
	void step_wrapper(Slots::Slots& slots, Slots::SlotsDevice& slots_d,
		const Background::Background& bkg, 
		const Background::BackgroundDevice& bkg_d, 
		const Options::Options& opts, const Options::OptionsDevice* opts_d,
		const OpenADAS::IonizRecombWarningsDevice& warnings_d, 
		int& max_time_ctr)
	{

#ifdef USE_CUDA
		if (opts.use_gpu_int() > 0) 
		{
			// Defined in impurity_transport.cu
			step_gpu(slots_d, opts.imp_time_step(), opts_d, warnings_d.max_time);
			return;
		}
#endif

		// Defined above
		step_cpu(slots, bkg, opts, opts.imp_time_step(), max_time_ctr);
	}


	/**
	* Updates particle charge states from ionization and recombination
	* probabilities, choosing the CPU or GPU implementation.
	*
	* The time step is the imp_time_step option. If GPU acceleration is on
	* (use_gpu), OpenADAS::ioniz_recomb_gpu is called. Otherwise
	* ioniz_recomb_cpu is called. The GPU branch only exists in CUDA builds.
	*
	* @param slots Host-side particle slots (used on the CPU path).
	* @param slots_d Device-side particle slots (used on the GPU path).
	* @param bkg Host-side background plasma (used on the CPU path).
	* @param bkg_d Device-side background plasma (used on the GPU path).
	* @param oa_ioniz Host-side OpenADAS ionization rate coefficients.
	* @param oa_ioniz_d Device-side OpenADAS ionization rate coefficients.
	* @param oa_recomb Host-side OpenADAS recombination rate coefficients.
	* @param oa_recomb_d Device-side OpenADAS recombination rate
	*                    coefficients.
	* @param opts Simulation options, which provide the time step.
	* @param ioniz_warnings Counter of particles whose ionization
	*                       probability exceeded 1.
	* @param recomb_warnings Counter of particles whose recombination
	*                        probability exceeded 1.
	* @param warnings_d Device-side warning counters (used on the GPU path,
	*                   where ioniz_warnings and recomb_warnings are not
	*                   updated).
	* @param rngs Host-side random number generators, one per OpenMP
	*             thread (used on the CPU path).
	* @param rngs_d Device-side random number generators (used on the GPU
	*               path).
	*/
	void ioniz_recomb_wrapper(Slots::Slots& slots, Slots::SlotsDevice& slots_d,
		const Background::Background& bkg, 
		const Background::BackgroundDevice& bkg_d, 
		const OpenADAS::OpenADAS& oa_ioniz,
		const OpenADAS::OpenADASDevice& oa_ioniz_d,
		const OpenADAS::OpenADAS& oa_recomb,
		const OpenADAS::OpenADASDevice& oa_recomb_d,
		const Options::Options& opts,
		int& ioniz_warnings, int& recomb_warnings,  
		const OpenADAS::IonizRecombWarningsDevice& warnings_d,
		std::vector<pcg32>& rngs, pcg32* rngs_d)
	{

#ifdef USE_CUDA
		if (opts.use_gpu_int() > 0) 
		{
			// Defined in openadas.cu
			OpenADAS::ioniz_recomb_gpu(slots_d, bkg_d, oa_ioniz_d, oa_recomb_d, 
				opts.imp_time_step(), warnings_d, rngs_d);
			return;
		}
#endif

		// Defined above
		ioniz_recomb_cpu(slots, bkg, oa_ioniz, oa_recomb, opts, 
			opts.imp_time_step(), ioniz_warnings, recomb_warnings, rngs);
	}


	/**
	* Applies Coulomb collisions with the background ions and electrons,
	* choosing the CPU or GPU implementation.
	*
	* If GPU acceleration is on (use_gpu), Collisions::collision_gpu is
	* called once for the ions and once for the electrons. The electron call
	* is skipped for the friction_force test case, which only considers ion
	* collisions. Otherwise collision_cpu is called, which makes the same
	* ion/electron distinction. The GPU branch only exists in CUDA builds.
	*
	* @param slots Host-side particle slots (used on the CPU path).
	* @param slots_d Device-side particle slots (used on the GPU path).
	* @param bkg Host-side background plasma (used on the CPU path).
	* @param bkg_d Device-side background plasma (used on the GPU path).
	* @param opts Simulation options.
	* @param imp_stats Host-side impurity statistics (used on the CPU
	*                  path).
	* @param imp_stats_d Device-side impurity statistics (used on the GPU
	*                    path).
	* @param dt Time step (s).
	* @param rngs Host-side random number generators, one per OpenMP
	*             thread (used on the CPU path).
	* @param rngs_d Device-side random number generators (used on the GPU
	*               path).
	* @param opts_d Device-side simulation options (used on the GPU path).
	*/
	void collision_wrapper(Slots::Slots& slots, 
		Slots::SlotsDevice& slots_d, const Background::Background& bkg, 
		const Background::BackgroundDevice& bkg_d, const Options::Options& opts,
		Impurity::Statistics& imp_stats, 
		ImpurityStats::StatisticsDevice& imp_stats_d, 
		const double dt, std::vector<pcg32>& rngs, pcg32* rngs_d, 
		const Options::OptionsDevice* opts_d)
	{

#ifdef USE_CUDA
		if (opts.use_gpu_int() > 0) 
		{
			// Defined in cuda/collision.cu. First call is for ions (the false)
			Collisions::collision_gpu(slots_d, bkg_d, false, dt, imp_stats_d,
				rngs_d, opts.gkyl_elec_mass_amu(), opts.gkyl_ion_mass_amu(), 
				opts_d);

			// friction_force test case only considers ion collisions to compare
			// against expected flow
			if (opts.test_opt_int() != 5)
			{
				Collisions::collision_gpu(slots_d, bkg_d, true, dt, imp_stats_d,
					rngs_d, opts.gkyl_elec_mass_amu(), opts.gkyl_ion_mass_amu(),
					opts_d);
			}
			return;
		}
#endif

		// The elec/ion distinction is done within
		collision_cpu(slots, bkg, dt, opts, imp_stats, rngs);
	}

	
	/**
	* Prints a warning for each warning counter that is nonzero.
	*
	* Warnings are printed for particles that exceeded the maximum allowed
	* time (imp_max_time), and for particles whose ionization or
	* recombination probability exceeded 1, which suggests the time step is
	* too large.
	*
	* @param opts Simulation options, used for imp_max_time in the warning
	*             message.
	* @param max_time_ctr Number of particles that exceeded the maximum
	*                     allowed time.
	* @param ioniz_warnings Number of particles whose ionization probability
	*                       exceeded 1.
	* @param recomb_warnings Number of particles whose recombination
	*                        probability exceeded 1.
	*/
	void print_counter_warnings(const Options::Options& opts, 
		int max_time_ctr, int ioniz_warnings, int recomb_warnings)
	{
		// Max time
		if (max_time_ctr)
			std::cout << "Warning: " << max_time_ctr << " particles exceeded "
				<< "maximum alotted time of " << opts.imp_max_time() << " s.\n";

		// Ionization warning
		if (ioniz_warnings)
			std::cout << "Warning: Ionization probability was > 1.0 for " 
				<< ioniz_warnings << " particles. Consider using a smaller "
				<< "time step.\n";

		// Recombination warning
		if (recomb_warnings)
			std::cout << "Warning: Recombination probability was > 1.0 for " 
				<< recomb_warnings << " particles. Consider using a smaller "
				<< "time step.\n";
	}


	/**
	* Runs the main particle-following loop until every particle has been
	* launched and followed until it dies.
	*
	* The imp_num particles are split evenly between the MPI ranks, with the
	* remainder going to the last rank. When running on GPUs, each rank's
	* share is split again between its GPUs, with the remainder going to
	* device 0. The slots are then filled, and the starting velocities are
	* pushed back half a time step as the Boris algorithm requires. Each
	* iteration of the loop then:
	*
	* 1. Updates charge states from ionization/recombination (if
	*    imp_iz_recomb is on).
	* 2. Applies collisions (if imp_collisions is on).
	* 3. Updates velocities with the Boris algorithm.
	* 4. Advances particle positions.
	* 5. Checks the time and spatial bounds.
	* 6. Replaces dead particles with new ones.
	* 7. Finds the containing cell of each particle.
	* 8. Records statistics.
	*
	* Each step is timed with a scoped timer. Progress is printed by rank 0
	* and GPU device 0 only.
	*
	* @param slots Host-side particle slots.
	* @param slots_d Device-side particle slots (GPU only).
	* @param bkg Host-side background plasma.
	* @param bkg_d Device-side background plasma (GPU only).
	* @param imp_stats Host-side impurity statistics.
	* @param imp_stats_d Device-side impurity statistics (GPU only).
	* @param oa_ioniz Host-side OpenADAS ionization rate coefficients.
	* @param oa_ioniz_d Device-side OpenADAS ionization rate coefficients
	*                   (GPU only).
	* @param oa_recomb Host-side OpenADAS recombination rate coefficients.
	* @param oa_recomb_d Device-side OpenADAS recombination rate
	*                    coefficients (GPU only).
	* @param opts Simulation options.
	* @param opts_d Device-side simulation options (GPU only).
	* @param timer Timer that accumulates the time spent in each step.
	* @param ioniz_warnings Counter of particles whose ionization
	*                       probability exceeded 1.
	* @param recomb_warnings Counter of particles whose recombination
	*                        probability exceeded 1.
	* @param warnings_d Device-side warning counters (GPU only). On the GPU
	*                   path, ioniz_warnings and recomb_warnings are not
	*                   updated here; the counts are read back at the end of
	*                   the run.
	* @param rngs Host-side random number generators, one per OpenMP
	*             thread (CPU only).
	* @param rngs_d Device-side random number generators (GPU only).
	* @param max_time_ctr Counter of particles that exceeded the maximum
	*                     allowed time.
	*/
	void main_loop(Slots::Slots& slots, Slots::SlotsDevice& slots_d,
		const Background::Background& bkg, 
		const Background::BackgroundDevice& bkg_d, 
		Impurity::Statistics& imp_stats, 
		ImpurityStats::StatisticsDevice& imp_stats_d,
		const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADASDevice& oa_ioniz_d, 
		const OpenADAS::OpenADAS& oa_recomb, 
		const OpenADAS::OpenADASDevice& oa_recomb_d, 
		Options::Options& opts, Options::OptionsDevice* opts_d, 
		Timer::Timer& timer, int& ioniz_warnings, int& recomb_warnings, 
		const OpenADAS::IonizRecombWarningsDevice& warnings_d,
		std::vector<pcg32>& rngs, pcg32* rngs_d, int& max_time_ctr)
	{

		// Rank and number of processes
		int rank {};
		int nprocs {};
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

		// Tracker for number of remaining particles to follow. We divide them
		// evenly between the number of MPI tasks/GPUs
		int parts_per_rank {opts.imp_num() / nprocs};

		// Give last process the remainder.
		if (rank == nprocs - 1) parts_per_rank += opts.imp_num() % nprocs;

		// Assign to remaining particles
		int rem_parts {parts_per_rank};
		
		// GPU device number
		int device {};

#ifdef USE_CUDA

		if (opts.use_gpu_int() > 0)
		{
			// Get device number, used below so only one GPU prints updates
			cudaGetDevice(&device);

			// Now subdivide among GPUs on this rank
			int num_gpus {};
			cudaGetDeviceCount(&num_gpus);

			int parts_per_gpu {rem_parts / num_gpus};
			int extra {rem_parts % num_gpus};

			// Evenly split between GPUs, assigning remainder to first one. I 
			// read somewhere that the first GPU is often the least loaded, but 
			// end of the day not a huge deal where the remainder goes since 
			// the number of GPUs isn't ever that large.
			if (slots_d.device_id == 0) 
				rem_parts = parts_per_gpu + extra;
			else rem_parts = parts_per_gpu;
		}
#endif

		// Struct to keep track of progress for user feedback.
		ParticleProgress prog(rem_parts, opts.print_interval());

		// Initial fill of slots with particles. alive_slots is used in
		// the progress print out in the loop
		int alive_slots {};
		fill_slots_wrapper(slots, slots_d, rem_parts, bkg, opts, alive_slots, 
			rngs, rngs_d, opts_d, bkg_d);

		// If for whatever reason the particle is started outside the grid,
		// check bounds to make sure it doesn't cause any chaos.
		check_bounds_wrapper(slots, slots_d, bkg, bkg_d, opts, rngs, rngs_d);

		// Find starting grid index
		find_containing_cell_wrapper(slots, slots_d, bkg, bkg_d, opts);

		// Boris algorithm: The velocity stored in the Slots objects
		// will actually be the velocity at a half timestep earlier, i.e.,
		// at t - dt/2. So we still need to push the particle velocity
		// back by half a time step.
		boris_wrapper(slots, slots_d, bkg, bkg_d, opts, 
			-opts.imp_time_step() / 2.0);

		// Record starting position in statistics arrays
		record_stats_wrapper(imp_stats, imp_stats_d, slots, slots_d, opts, 
			opts.imp_time_step());

		// Begin loop
		bool all_dead {false};
		while (!all_dead)
		{

			// We create a scope for each step so that we can use a scoped
			// timer to profile the time spent in each step of the loop.
			// A scoped timer is a bit safer because it automatically stops
			// when it goes out of scope.

#ifdef DEBUG
			std::cout << "Wrapper called: Ionization/recombination...\n";
#endif

			// Check for ionization/recombination
			if (opts.imp_iz_recomb_int() > 0)
			{
				Timer::ScopedTimer t(timer.acc(Timer::Section::IonRec));
				ioniz_recomb_wrapper(slots, slots_d, bkg, bkg_d, oa_ioniz,
					oa_ioniz_d, oa_recomb, oa_recomb_d, opts, ioniz_warnings,
					recomb_warnings, warnings_d, rngs, rngs_d);
			}

			// Variance reduction
			// To-do

#ifdef DEBUG
			std::cout << "Wrapper called: Collisions...\n";
#endif

			// Collision update
			if (opts.imp_collisions_int() > 0)
			{
				Timer::ScopedTimer t(timer.acc(Timer::Section::Coll));
				collision_wrapper(slots, slots_d, bkg, bkg_d, opts,
					imp_stats, imp_stats_d, opts.imp_time_step(), rngs, rngs_d,
					opts_d);
			}

#ifdef DEBUG
			std::cout << "Wrapper called: Boris...\n";
#endif

			// Update velocity (Boris).
			{
				Timer::ScopedTimer t(timer.acc(Timer::Section::Boris));
				boris_wrapper(slots, slots_d, bkg, bkg_d, opts, 
					opts.imp_time_step());
			}

#ifdef DEBUG
			std::cout << "Wrapper called: Step...\n";
#endif

			// Perform particle step
			{
				Timer::ScopedTimer t(timer.acc(Timer::Section::Step));
				step_wrapper(slots, slots_d, bkg, bkg_d, opts, opts_d, 
					warnings_d, max_time_ctr);
			}

#ifdef DEBUG
			std::cout << "Wrapper called: Bounds...\n";
#endif

			// Bounds checking
			{
				Timer::ScopedTimer t(timer.acc(Timer::Section::Bounds));
				check_bounds_wrapper(slots, slots_d, bkg, bkg_d, opts, rngs, 
				rngs_d);
			}

#ifdef DEBUG
			std::cout << "Wrapper called: Fill slots...\n";
#endif

			// Replace dead particles
			{
				Timer::ScopedTimer t(timer.acc(Timer::Section::FillSlots));
				fill_slots_wrapper(slots, slots_d, rem_parts, bkg, opts, 
					alive_slots, rngs, rngs_d, opts_d, bkg_d);
			}

#ifdef DEBUG
			std::cout << "Wrapper called: Find cell...\n";
#endif

			// Update particle indices
			{
				Timer::ScopedTimer t(timer.acc(Timer::Section::FindCell));
				find_containing_cell_wrapper(slots, slots_d, bkg, bkg_d, opts);
			}

#ifdef DEBUG
			std::cout << "Wrapper called: Record statistics...\n";
#endif

			// Record statistics
			{
				Timer::ScopedTimer t(timer.acc(Timer::Section::Record));
				record_stats_wrapper(imp_stats, imp_stats_d, slots, slots_d, 
					opts, opts.imp_time_step());
			}

			// User feedback just for rank 0, no MPI communication. The total
			// number of remaining particles is rem_parts plus the number
			// actively being followed (alive_slots).
			if (rank == 0 && device == 0) 
				prog.update(rem_parts + alive_slots);

#ifdef DEBUG
			std::cout << "Wrapper called: All dead...\n";
#endif

			// Check if all the particles in slots are dead. If this happens
			// after fill_slots, it means there were no more alive particles
			// to swap in and all the remaining ones are dead. So we're done.
			// For GPUs, this has a GPU->CPU memory copy which means it has
			// an implicit synchronize in it. 
			{
				Timer::ScopedTimer t(timer.acc(Timer::Section::AllDead));
				all_dead = all_dead_wrapper(slots, slots_d, opts, rem_parts);
			}

			// Intermittent save
			// To-do

		}  // while (!all_dead)
	}  // main_loop

#ifdef USE_CUDA

	/**
	* Runs main_loop on a single GPU from its own host thread.
	*
	* follow_impurities starts one of these per GPU so the GPUs run
	* concurrently. The thread first selects the GPU given by
	* slots_d.device_id, then calls main_loop with that GPU's device-side
	* data. An empty host-side RNG vector is passed because it is not used
	* on the GPU path. Only defined in CUDA builds.
	*
	* @param slots Host-side particle slots.
	* @param slots_d Device-side particle slots for this GPU.
	* @param bkg Host-side background plasma.
	* @param bkg_d Device-side background plasma for this GPU.
	* @param imp_stats Host-side impurity statistics.
	* @param imp_stats_d Device-side impurity statistics for this GPU.
	* @param oa_ioniz Host-side OpenADAS ionization rate coefficients.
	* @param oa_ioniz_d Device-side OpenADAS ionization rate coefficients
	*                   for this GPU.
	* @param oa_recomb Host-side OpenADAS recombination rate coefficients.
	* @param oa_recomb_d Device-side OpenADAS recombination rate
	*                    coefficients for this GPU.
	* @param opts Simulation options.
	* @param timer Timer for this GPU.
	* @param ioniz_warnings Counter of particles whose ionization
	*                       probability exceeded 1.
	* @param recomb_warnings Counter of particles whose recombination
	*                        probability exceeded 1.
	* @param warnings_d Device-side ionization/recombination warning
	*                   counters for this GPU.
	* @param rngs_d Device-side random number generators for this GPU.
	* @param opts_d Device-side simulation options for this GPU.
	* @param max_time_ctr Counter of particles that exceeded the maximum
	*                     allowed time.
	*/
	void gpu_worker(Slots::Slots& slots, Slots::SlotsDevice& slots_d,
		const Background::Background& bkg, 
		const Background::BackgroundDevice& bkg_d,
		Impurity::Statistics& imp_stats, 
		ImpurityStats::StatisticsDevice& imp_stats_d,
		const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADASDevice& oa_ioniz_d, 
		const OpenADAS::OpenADAS& oa_recomb, 
		const OpenADAS::OpenADASDevice& oa_recomb_d, 
		Options::Options& opts, Timer::Timer& timer, int& ioniz_warnings, 
		int& recomb_warnings, 
		const OpenADAS::IonizRecombWarningsDevice& warnings_d,
		pcg32* rngs_d, Options::OptionsDevice* opts_d,
		int& max_time_ctr)
	{
		// Not used, so just dummying
		std::vector<pcg32> rngs;

		cudaSetDevice(slots_d.device_id);   // THIS THREAD USES THIS GPU
		main_loop(slots, slots_d, bkg, bkg_d, imp_stats, imp_stats_d, oa_ioniz, 
			oa_ioniz_d, oa_recomb, oa_recomb_d, opts, opts_d, timer, 
			ioniz_warnings, recomb_warnings, warnings_d, rngs, rngs_d, 
			max_time_ctr);
	}

#endif

	/**
	* Follows every impurity particle through the background plasma and
	* returns the accumulated statistics.
	*
	* This is the top-level driver for impurity transport. It sets up the
	* impurity statistics arrays, loads the OpenADAS rate coefficients (one
	* MPI rank at a time, to avoid every rank reading the same files at
	* once) and allocates the particle slots, with min(slot_cap, imp_num)
	* slots. All slots start dead so the first fill launches particles into
	* them.
	*
	* On GPUs, the data is copied to each GPU, a random number generator is
	* set up for each GPU thread, and gpu_worker runs main_loop on every GPU
	* concurrently. The statistics and timers from each GPU are then combined
	* and the device memory is freed. On CPUs, one random number generator
	* is created per OpenMP thread and main_loop is called directly. If GPU
	* acceleration was requested in a build without CUDA, it is turned off.
	*
	* If the seed option is negative, random number generators are seeded
	* from the clock so each run is different. Otherwise the seed option is
	* used. In both cases the MPI rank is added to the seed so each rank
	* gets a different sequence.
	*
	* @param bkg Background plasma to follow the particles through.
	* @param opts Simulation options. use_gpu may be turned off if Flan was
	*             built without GPU support.
	* @param timer Timer that accumulates the time spent in each step.
	* @return Impurity statistics accumulated from every particle followed
	*         on this rank.
	*/
	Impurity::Statistics follow_impurities(Background::Background& bkg, 
		Options::Options& opts, Timer::Timer& timer)
	{
		
#ifdef DEBUG
			std::cout << "Debug mode enabled. Synchronizations after each" 
				<< "wrapper call and print statements are on.\n";
#endif

		// Some warnings to include as I still implement things
		std::cout << "Warning! Need to still implement lcfs_x to divide"
			<< " core/SOL BCs\n";

		// Rank and number of processes
		int rank {};
		int nprocs {};
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

		// Initialize particle statistics vectors, all contained within a
		// Statistics object. 
		Impurity::Statistics imp_stats {bkg.get_dim1(), bkg.get_dim2(), 
			bkg.get_dim3(), bkg.get_dim4()};

		// Load OpenADAS data needed for ionization/recombination rates if
		// ionization/recombination is on. Believe it or not, this is the
		// preferred way to do this in C++.
		// We only want one MPI process to do this at a time since it involves
		// reading the ADAS files, and having different processes all reading
		// the same file could cause slow downs. Really minor impact though.
		OpenADAS::OpenADAS oa_ioniz {};
		OpenADAS::OpenADAS oa_recomb {};
		for (int r {}; r < nprocs; r++)
		{
			if (rank == r)
			{
				oa_ioniz = opts.imp_iz_recomb_int() == 1 
					? OpenADAS::OpenADAS(opts.openadas_root(), 
					opts.openadas_year(), opts.imp_atom_num(), "scd") 
					: OpenADAS::OpenADAS();
				oa_recomb = opts.imp_iz_recomb_int() == 1 
					? OpenADAS::OpenADAS(opts.openadas_root(), 
					opts.openadas_year(), opts.imp_atom_num(), "acd") 
					: OpenADAS::OpenADAS();
			}

			// Everyone waits before next rank runs
			MPI_Barrier(MPI_COMM_WORLD);
		}

		// Used to keep track of when the ionization/recombination probability
		// is greater than 1 (from too large a time step).
		int ioniz_warnings {};
		int recomb_warnings {};

		// Counter to see how many particles exceed the maximum alotted time
		// (in plasma time)
		int max_time_ctr {};

		// Slot capacity, 2^20 seems reasonable but can be adjusted. The larger
		// the better (probably). The motivation for this is to alleviate 
		// memory pressure, since you could theoretically want to follow
		// millions or billions of particles, and allocating them all at once
		// is wasteful and could cause you to run out of memory.
		//constexpr int slot_cap = 1048576;  
		int nslots {std::min(opts.slot_cap(), opts.imp_num())};
		std::cout << "Slot Capacity = " << nslots << '\n';
		Slots::Slots slots {nslots};

		// All particles assumed to have same mass and atomic number.
		slots.set_mass(opts.imp_mass_amu() * Constants::amu_to_kg);
		slots.set_Z(opts.imp_atom_num());

		// Set all particles to dead so that we can pass through the slot
		// filling logic that will assign initial conditions
		for (int i {}; i < slots.N(); ++i)
		{
			slots.set_state(i, 1);
		}

#ifdef USE_CUDA

		// Get number of GPUs
		int num_gpus {};
		if (opts.use_gpu_int() > 0)
		{
			cudaGetDeviceCount(&num_gpus);
			std::cout << "Using " << num_gpus << " GPU(s)\n";
		}
		else
		{
			std::cout << "GPU acceleration is off. Turn on with "
				<< "use_gpu = \"cuda\"\n";
		}

		if (opts.use_gpu_int() > 0)
		{

#ifdef DEBUG
			std::cout << "Copying data to GPU(s)...\n";
#endif

			// Create device-side structs on each GPU if using GPUs
			std::vector<Slots::SlotsDevice> gpu_slots;
			std::vector<ImpurityStats::StatisticsDevice> gpu_stats;
			std::vector<Background::BackgroundDevice> gpu_bkgs;
			std::vector<OpenADAS::OpenADASDevice> gpu_oa_izs;
			std::vector<OpenADAS::OpenADASDevice> gpu_oa_rcs;
			std::vector<Options::OptionsDevice*> gpu_opts;
			std::vector<pcg32*> gpu_rngs;
			std::vector<OpenADAS::IonizRecombWarningsDevice> gpu_warnings;
			std::vector<Timer::Timer> gpu_timers;
			for (int dev = 0; dev < num_gpus; dev++) 
			{
				gpu_slots.push_back(slots.to_device(dev));
				gpu_stats.push_back(imp_stats.to_device(dev));
				gpu_bkgs.push_back(bkg.to_device(dev));
				gpu_oa_izs.push_back(oa_ioniz.to_device(dev));
				gpu_oa_rcs.push_back(oa_recomb.to_device(dev));
				gpu_warnings.push_back(OpenADAS::alloc_warnings(dev));

				// Note that we are storing OptionsDevice* here. It's treated
				// differently since there aren't arrays to be copied, so 
				// instead of copying a struct of multiple device-side pointers
				// each kernel call, we copy a single pointer to a device-side
				// struct. I think this optimization is negligible, but
				// whatever had to pick something.
				gpu_opts.push_back(opts.to_device(dev));

				// Initialize RNGs for each GPU thread. Uses the just
				// allocated slots on the GPU. gpu_rngs is a vector of PCG32
				// arrays. 
				// If seed is negative (default), initialize the base seed from the
				// clock so each run is different. MPI ranks are distinguished 
				// through the seed, threads are distinguished through the 
				// PCG stream, and GPU are distinguished by the number (dev).
				uint64_t seed {};
				if (opts.seed() < 0)
				{
					// Nanoseconds since epoch, mixed with rank and thread ID
					seed = static_cast<uint64_t>(
						std::chrono::high_resolution_clock::now().time_since_epoch()
						.count()) + rank * num_gpus + dev;
				}
				else
				{
						seed = static_cast<uint64_t>(opts.seed()) + rank 
							* num_gpus + dev;
				}

				gpu_rngs.push_back(init_rngs_cuda(gpu_slots[dev], seed));

				// Each GPU gets it own timer
				gpu_timers.push_back(Timer::Timer());
			}

#ifdef DEBUG
			std::cout << "Starting main loop (GPU)...\n";
#endif

			// Spawn threads to launch main_loop on each available GPU
			std::vector<std::thread> threads;
			for (int dev = 0; dev < num_gpus; dev++) 
			{

				// Start a GPU worker on each thread. std::ref is used here 
				// because std::thread copies by default, which isn't necessary
				// or desirable for these larger objects.
				threads.emplace_back(gpu_worker, std::ref(slots), 
					std::ref(gpu_slots[dev]), std::ref(bkg), 
					std::ref(gpu_bkgs[dev]), std::ref(imp_stats), 
					std::ref(gpu_stats[dev]), 
					std::ref(oa_ioniz), std::ref(gpu_oa_izs[dev]), 
					std::ref(oa_recomb), std::ref(gpu_oa_rcs[dev]), 
					std::ref(opts), std::ref(gpu_timers[dev]), 
					std::ref(ioniz_warnings), std::ref(recomb_warnings), 
					std::ref(gpu_warnings[dev]),
					std::ref(gpu_rngs[dev]), std::ref(gpu_opts[dev]), 
					std::ref(max_time_ctr));
			}

			// Wait for all threads to finish
			for (auto& t : threads) {
				t.join();
			}

			for (int dev = 0; dev < num_gpus; dev++) 
			{
				// Reduce GPU stats
				std::cout << "Reducing stats...\n";
				imp_stats.add_stats_device(gpu_stats[dev], dev);

				// Combine timers into initial one
				timer.merge(gpu_timers[dev]);

				// Free memory of device-side structs
				Slots::free_slots(gpu_slots[dev], dev);
				Impurity::free_stats(gpu_stats[dev], dev);
				Background::free_bkg(gpu_bkgs[dev], dev);
				OpenADAS::free_oa(gpu_oa_izs[dev], dev);
				OpenADAS::free_oa(gpu_oa_rcs[dev], dev);
				Options::free_opts(gpu_opts[dev], dev);
				free_rngs(gpu_rngs[dev], dev);

				// Add up the warnings from each GPU and free the counters.
				// Done here, after the threads have joined, so that only one
				// thread ever writes to ioniz_warnings/recomb_warnings.
				OpenADAS::add_warnings_to_host(gpu_warnings[dev], 
					ioniz_warnings, recomb_warnings, max_time_ctr);
				OpenADAS::free_warnings(gpu_warnings[dev]);
			}

			print_counter_warnings(opts, max_time_ctr, ioniz_warnings, 
				recomb_warnings);
			return imp_stats;
		}

#endif

		// Dummy device-side structs. They're not used in CPU-only simulations.
		Slots::SlotsDevice slots_d {};
		ImpurityStats::StatisticsDevice imp_stats_d {};
		Background::BackgroundDevice bkg_d {};
		OpenADAS::OpenADASDevice oa_ioniz_d {};
		OpenADAS::OpenADASDevice oa_recomb_d {};
		OpenADAS::IonizRecombWarningsDevice warnings_d {};
		Options::OptionsDevice* opts_d {};
		pcg32* rngs_d {};

		// Vector of PCG32 RNGs, one per OpenMP thread per MPI rank. We want
		// to pass these through all the various functions that need random
		// numbers. It's important that each thread has its own RNG to avoid
		// overlapping sequences. 
		std::vector<pcg32> rngs;
		int nthreads = omp_get_max_threads();
		rngs.reserve(nthreads);

		for (int tid = 0; tid < nthreads; ++tid)
		{
			uint64_t seed {};
			uint64_t stream {static_cast<uint64_t>(tid)}; // unique per thread
			
			// If seed is negative (default), initialize the base seed from the
			// clock so each run is different. MPI ranks are distinguished 
			// through the seed, while threads are distinguished through the 
			// PCG stream.
			if (opts.seed() < 0)
			{
				// Nanoseconds since epoch, mixed with rank and thread ID
				seed = static_cast<uint64_t>(
					std::chrono::high_resolution_clock::now().time_since_epoch()
					.count()) + rank;
			}
			else
			{
					seed = static_cast<uint64_t>(opts.seed()) + rank;
			}
		
			rngs.push_back(pcg32(seed, stream));
		}

#ifdef DEBUG
			std::cout << "Starting main loop (CPU)...\n";
#endif

		main_loop(slots, slots_d, bkg, bkg_d, imp_stats, imp_stats_d, oa_ioniz, 
			oa_ioniz_d, oa_recomb, oa_recomb_d, opts, opts_d, timer, 
			ioniz_warnings, recomb_warnings, warnings_d, rngs, rngs_d, 
			max_time_ctr);

		print_counter_warnings(opts, max_time_ctr, ioniz_warnings, 
			recomb_warnings);

		return imp_stats;
	}  // follow_impurities
}  // namespace ImpurityTransport
