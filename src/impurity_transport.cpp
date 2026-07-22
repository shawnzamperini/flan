#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <thread>
#include <tuple>
#include <vector>

#include "background.h"
#include "boris.h"
#include "collisions.h"
#include "constants.h"
#include "flan_types.h"
#include "impurity.h"
#include "impurity_stats.h"
#include "impurity_transport.h"
#include "openadas.h"
#include "options.h"
#include "random.h"
#include "slots.h"
#include "timer.h"
#include "utilities.h"
#include "variance_reduction.h"
#include "vectors.h"

#ifdef USE_CUDA
#include <cuda_runtime.h>
#include "impurity_stats.cuh"
#include "impurity_transport.cuh"
#include "slots.cuh"
#endif


namespace ImpurityTransport
{

	template <typename T>
	int get_nearest_index_cpu(const std::vector<T>& vec, const T value)
	{
		auto lower = std::lower_bound(vec.begin(), vec.end(), value);
    
		if (lower == vec.begin()) 
		{
			return 0;
		}
		if (lower == vec.end()) 
		{
			return vec.size() - 1;
		}
		
		auto prev = lower - 1;
		if (std::fabs(*lower - value) < std::fabs(*prev - value)) 
		{
			return lower - vec.begin();
		} 
		else 
		{
			return prev - vec.begin();
		}
	}

	template <typename T>
	int get_nearest_cell_index_cpu(const std::vector<T>& grid_edges, 
		const T value)
	{
		// Get the index of the first value in grid_edges that is larger
		// than value.
		auto lower = std::lower_bound(grid_edges.begin(), grid_edges.end(), 
			value);

		// Realize that one minus the index represented by lower is the value
		// we're after in the vectors with values at the cell centers.
		//  ____________
		//  |_0_|_1_|_2_|  <-- cell center indices
		//  0   1   2   3  <-- grid_edges indices
		//          ^
		//        lower
		//
		// In this example, we want 1 returned, so we return 2 - 1 = 1. If
		// value is outside the range, return the index of the respective end
		// of the vector.
		int index = std::distance(grid_edges.begin(), lower);

		// If less than everything, return the first cell.
		if (lower == grid_edges.begin()) return 0;

		// If larger than everything, return the last cell. Note end() is the
		// iterator that points past the last element of a vector, so to return
		// the cell we need to subtract by 2. 
		else if (lower == grid_edges.end()) 
		{
			return index - 2;
		}

		//else return lower - grid_edges.begin() - 1;
		else return index - 1;
	}

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
			int tidx {get_nearest_index_cpu(times, slots.t()[i])};

			// Find nearest cell using the grid (unlike time, which uses
			// cell "centers"). 
			int xidx {get_nearest_cell_index_cpu(grid_x, slots.x()[i])};
			int yidx {get_nearest_cell_index_cpu(grid_y, slots.y()[i])};
			int zidx {get_nearest_cell_index_cpu(grid_z, slots.z()[i])};

			slots.set_tidx(i, tidx);
			slots.set_xidx(i, xidx);
			slots.set_yidx(i, yidx);
			slots.set_zidx(i, zidx);
		}
	}

	// Absorbing boundary condition. Check is value is outside the bound.
	// max_bound = true indicates that bound is a maximum, and that we are
	// checking if value >= bound. Likewise, max_bound = false checks if
	// value =< bound. buffer extends the boundary so that it triggers sooner.
	// This is because sometimes the plasma can behave oddly right at the
	// boundaries, so this can circumvent it. This function is written this way 
	// to be SIMD/vectorizable friendly for the compiler (no if's).
	inline int absorbing_bc(double a, double bound, bool max_bound,
		double buffer = 0.0)
	{
		// Compute both comparisons
		const bool ge {(a + buffer) >= bound};  // for max-bound
		const bool le {(a - buffer) <= bound};  // for min-bound

		// Select the correct one using max_bound as a mask
		const bool dead {max_bound ? ge : le};

		// Returns 0 if alive, 1 if dead
		return static_cast<int>(dead);
	}

	// Periodic boundary condition. This type of BC will never kill a particle,
	// just loop it around to the other side. This function returns the new
	// looped value, or just the same value if a BC was not encountered.
	inline double periodic_bc(double a, double amin, double amax, 
		double buffer = 0.0)
	{
		const double L = amax - amin;

		const double gt = ((a + buffer) > amax);   // wrap high, subtract L
		const double lt = ((a - buffer) < amin);   // wrap low, add L

		// Return new (or same) position
		return a - gt * L + lt * L;
	}

	// Core boundary condition. The coordinate being checked against is a. b
	// and c are the other two coordinates. If a core condition is met, b is
	// moved to a random value between b_min/b_max, likewise for c. This is
	// essentially like entering the core and then popping out at a random
	// location somewhere else along the core boundary. A residence time
	// could easily be added to this if it was useful.
	inline void core_bc(double& a, double& b, double& c, double a_min,
		double a_max, double buffer, double b_min, double b_max,
		double c_min, double c_max)
	{
		// Branchless boundary detection
		const double hit_min = (a - buffer) < a_min;
		const double hit_max = (a + buffer) > a_max;

		// Combined mask: 1.0 if either boundary is hit
		const double mask = hit_min || hit_max;

		// Compute remap positions for min and max sides
		const double a_new_min = a_min + buffer;
		const double a_new_max = a_max - buffer;

		// Select correct remap position
		const double a_new =
			hit_min * a_new_min +
			hit_max * a_new_max;

		// Branchless update of the coordinate
		a = mask * a_new + (1.0 - mask) * a;

		// RNG only when needed (cannot be SIMD-friendly)
		if (mask)
		{
			b = Random::get(b_min, b_max);
			c = Random::get(c_min, c_max);
		}
	}


	void check_bounds_cpu(Slots::Slots& slots, 
		const Background::Background& bkg, const Options::Options& opts)
	{

		// The different bounds are
		//   0: Absorbing
		//   1: Periodic
		//   2: Core (only x,y,z)

		// Loop through each slot
		#pragma omp parallel for
		for (int i = 0; i < slots.N(); i++)
		{
			// Skip dead particles
			if (slots.state()[i] > 0) continue;

			// Used in all, just define up front
			int state {};
			double new_val {};

			// --------------------
			// Time boundary
			// --------------------

			// Maximum t: Absorbing boundary
			if (opts.tbound_type_int() == 0)
			{
				state = absorbing_bc(slots.t()[i], bkg.get_t_max(), true);
				slots.set_state(i, state);

				// If dead (>0) we can skip the rest of the checks.
				if (state) continue;
			}

			// Maximum t: Periodic boundary
			else if (opts.tbound_type_int() == 1)
			{
				// We check and overwrite each time here, even though hitting
				// a BC is relatively rare. Despite this, this is still a very
				// fast process since slots.t()[i] is loaded into L1, and 
				// periodic_bc doesn't touch memory loads so the L1 cache never
				// get evicted, so it's a very fast write.
				new_val = periodic_bc(slots.t()[i], bkg.get_t_min(), 
					bkg.get_t_max());
				slots.set_t(i, new_val);
			}

			// --------------------
			// Minimum x boundary
			// --------------------
		
			// Absorbing boundary
			if (opts.min_xbound_type_int() == 0)
			{
				state = absorbing_bc(slots.x()[i], bkg.get_x_min(), false, 
					opts.imp_xbound_buffer());
				slots.set_state(i, state);

				// If dead (>0) we can skip the rest of the checks.
				if (state) continue;
			}

			// Periodic boundary
			else if (opts.min_xbound_type_int() == 1)
			{
				new_val = periodic_bc(slots.x()[i], bkg.get_x_min(), 
					bkg.get_x_max(), opts.imp_xbound_buffer());
				slots.set_x(i, new_val);
			}

			// Core boundary
			else if (opts.min_xbound_type_int() == 2)
			{

				// Temporaries that may get updated in core_bc, which we will
				// write into slots after
				double xtmp {slots.x()[i]};
				double ytmp {slots.y()[i]};
				double ztmp {slots.z()[i]};
				core_bc(xtmp, ytmp, ztmp, bkg.get_x_min(), bkg.get_x_max(), 
					opts.imp_xbound_buffer(), bkg.get_y_min(), bkg.get_y_max(), 
					bkg.get_z_min(), bkg.get_z_max());

				// Write back
				slots.set_x(i, xtmp);
				slots.set_y(i, ytmp);
				slots.set_z(i, ztmp);
			}

			// --------------------
			// Maximum x boundary
			// --------------------

			// Absorbing boundary
			if (opts.max_xbound_type_int() == 0)
			{
				state = absorbing_bc(slots.x()[i], bkg.get_x_max(), true, 
					opts.imp_xbound_buffer());
				slots.set_state(i, state);

				// If dead (>0) we can skip the rest of the checks.
				if (state) continue;
			}

			// Periodic boundary
			else if (opts.max_xbound_type_int() == 1)
			{
				new_val = periodic_bc(slots.x()[i], bkg.get_x_min(), 
					bkg.get_x_max(), opts.imp_xbound_buffer());
				slots.set_x(i, new_val);
			}

			// Core boundary
			else if (opts.max_xbound_type_int() == 2)
			{

				// Temporaries that may get updated in core_bc, which we will
				// write into slots after
				double xtmp {slots.x()[i]};
				double ytmp {slots.y()[i]};
				double ztmp {slots.z()[i]};
				core_bc(xtmp, ytmp, ztmp, bkg.get_x_min(), bkg.get_x_max(), 
					opts.imp_xbound_buffer(), bkg.get_y_min(), bkg.get_y_max(), 
					bkg.get_z_min(), bkg.get_z_max());

				// Write back
				slots.set_x(i, xtmp);
				slots.set_y(i, ytmp);
				slots.set_z(i, ztmp);
			}


			// --------------------
			// y boundary
			// --------------------

			// Minimum y: Periodic y boundary
			if (slots.y()[i] < bkg.get_grid_y()[0])
			{
				slots.set_y(i, bkg.get_grid_y().back() + (slots.y()[i] 
					- bkg.get_grid_y()[0]));
			}
			else if (slots.y()[i] > bkg.get_grid_y().back())
			{
				slots.set_y(i, bkg.get_grid_y()[0] + (slots.y()[i]
					- bkg.get_grid_y().back()));
			}

			// --------------------
			// z boundary
			// --------------------

			// Absorbing z boundary in SOL, periodic in core
			// Problem: This assume x increases with distance from the core,
			// this is not always true! E.g., it could depend on how psi is
			// defined. 
			if (slots.x()[i] > opts.lcfs_x())
			{
				// Minimum/maximum z: Absorbing boundary in SOL
				if (slots.z()[i] < bkg.get_grid_z()[0] || 
					slots.z()[i] > bkg.get_grid_z().back()) 
					{
						//std::cout << "Absorbed: Min/max z\n";
						slots.set_state(i, 1);
						continue;
					}
			}
			else
			{
				// Minimum z: Periodic boundary
				if (slots.z()[i] < bkg.get_grid_z()[0])
				{
					slots.set_z(i, bkg.get_grid_z().back() + (slots.z()[i] 
						- bkg.get_grid_z()[0]));
					//std::cout << "Periodic: Min z\n";
				}
				
				// Maximum z: Periodic boundary
				else if (slots.z()[i] > bkg.get_grid_z().back())
				{
					slots.set_z(i, bkg.get_grid_z()[0] + (slots.z()[i] 
						- bkg.get_grid_z().back()));
					//std::cout << "Periodic: Max z\n";
				}
			}

		} // slot loop
	}  // check_bounds_cpu


	// Wrapper to choose CPU or GPU implementation for finding containing cell
	// for each particle in a slot
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


	// Wrapper to choose CPU or GPU implementation for recording particle
	// statistics in the underlying statistics arrays.
	void record_stats_wrapper(Impurity::Statistics& imp_stats, 
		ImpurityStats::StatisticsDevice& imp_stats_d, 
		const Slots::Slots& slots, const Slots::SlotsDevice& slots_d, 
		const Options::Options& opts, const double imp_time_step)
	{

#ifdef USE_CUDA
		if (opts.use_gpu_int() > 0)
		{
			// Defined in cuda/impurity_transport.cu
			ImpurityStats::record_stats_gpu(imp_stats_d, slots_d, imp_time_step);
			return;
		}
#endif

		// Defined in impurity_stats.cpp
		record_stats_cpu(imp_stats, slots, opts, imp_time_step);
	}


	void fill_slots_wrapper(Slots::Slots& slots, Slots::SlotsDevice& slots_d,
		int& rem_parts, Options::Options& opts)
	{

#ifdef USE_CUDA
		if (opts.use_gpu_int() > 0) 
		{
			// Defined in cuda/slots.cu
			Slots::fill_slots_gpu(slots_d, rem_parts);
			return;
		}
#endif

		// Defined in slots.cpp
		Slots::fill_slots_cpu(slots, rem_parts);
	}


	void boris_wrapper(Slots::Slots& slots, 
		Slots::SlotsDevice& slots_d, const Background::Background& bkg, 
		const Background::BackgroundDevice& bkg_d, const Options::Options& opts,
		const double dt)
	{

#ifdef USE_CUDA
		if (opts.use_gpu_int() > 0) 
		{
			// Boris::update_velocity_gpu
			return;
		}
#endif

		Boris::update_velocity_cpu(slots, slots_d, bkg, bkg_d, opts, dt);

	}

	// Check if all particles in slots are dead
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


	// Check if time or spatial bounds have been exceeded and handle 
	// appropriately
	void check_bounds_wrapper(Slots::Slots& slots, Slots::SlotsDevice& slots_d,
		const Background::Background& bkg, 
		const Background::BackgroundDevice& bkg_d, 
		const Options::Options& opts)
	{

#ifdef USE_CUDA
		if (opts.use_gpu_int() > 0) 
		{
			ImpurityTransport::check_bounds_gpu(slots_d, bkg_d, 
				opts.tbound_type_int(), opts.imp_xbound_buffer(), 
				opts.min_xbound_type_int(), opts.lcfs_x());
			return;
		}
#endif

		check_bounds_cpu(slots, bkg, opts);
	}

	
	void main_loop(Slots::Slots& slots, 
		Slots::SlotsDevice& slots_d,
		const Background::Background& bkg, 
		const Background::BackgroundDevice& bkg_d, 
		Impurity::Statistics& imp_stats, 
		ImpurityStats::StatisticsDevice& imp_stats_d,
		const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADAS& oa_recomb, Options::Options& opts, 
		Timer::Timer& timer)
	{

		// Rank and number of processes
		int rank {};
		int nprocs {};
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

		// Tracker for number of remaining particles to follow. We divide them
		// evenly between the number of MPI tasks/GPUs
		//int rem_parts {opts.imp_num()};
		int parts_per_rank {opts.imp_num() / nprocs};

		// Give last process the remainder.
		if (rank == nprocs - 1) parts_per_rank += opts.imp_num() % nprocs;

		// Assign to remaining particles
		int rem_parts {parts_per_rank};
		
#ifdef USE_CUDA

		if (opts.use_gpu_int() > 0)
		{
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

		// Initial fill of slots with particles. 
		std::cout << "pre-fill: rem_parts = " << rem_parts << '\n';
		fill_slots_wrapper(slots, slots_d, rem_parts, opts);
		std::cout << "post-fill: rem_parts = " << rem_parts << '\n';
		
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

			// Check for ionization/recombination
			// To-do

			// Variance reduction
			// To-do

			// Collision update

			// Update velocity (Boris)
			boris_wrapper(slots, slots_d, bkg, bkg_d, opts, 
				opts.imp_time_step());

			// Perform particle step
			// Temporary just to debug
			for (int i=0; i < slots.N(); ++i)
			{
				slots.set_t(i, slots.t()[i] + opts.imp_time_step());
			}

			// Bounds checking
			check_bounds_wrapper(slots, slots_d, bkg, bkg_d, opts);

			// Update particle indices
			find_containing_cell_wrapper(slots, slots_d, bkg, bkg_d, opts);

			// Record statistics
			record_stats_wrapper(imp_stats, imp_stats_d, slots, slots_d, opts, 
				opts.imp_time_step());

			// Replace dead particles
			std::cout << "pre-fill: rem_parts = " << rem_parts << '\n';
			fill_slots_wrapper(slots, slots_d, rem_parts, opts);
			std::cout << "post-fill: rem_parts = " << rem_parts << '\n';

			// User feedback

			// Check if all the particles in slots are dead. If this happens
			// after fill_slots, it means there were no more alive particles
			// to swap in and all the remaining ones are dead. So we're done.
			all_dead = all_dead_wrapper(slots, slots_d, opts, rem_parts);

			// Intermittent save
			// To-do

		}  // while (!all_dead)
	
	
	}

#ifdef USE_CUDA

	// Function to enable each thread to control each GPU independently to 
	// allow them to run concurrently. Used when multiple GPUs are available.
	void gpu_worker(Slots::Slots& slots, Slots::SlotsDevice& slots_d,
		const Background::Background& bkg, 
		const Background::BackgroundDevice& bkg_d,
		Impurity::Statistics& imp_stats, 
		ImpurityStats::StatisticsDevice imp_stats_d,
		const OpenADAS::OpenADAS& oa_ioniz, const OpenADAS::OpenADAS& oa_recomb, 
		Options::Options& opts, Timer::Timer& timer)
	{
		cudaSetDevice(slots_d.device_id);   // THIS THREAD USES THIS GPU
		main_loop(slots, slots_d, bkg, bkg_d, imp_stats, imp_stats_d, oa_ioniz, 
			oa_recomb, opts, timer);
	}

#endif

	Impurity::Statistics follow_impurities(Background::Background& bkg, 
		Options::Options& opts, Timer::Timer& timer)
	{
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

		// Slot capacity, 2^17 seems reasonable but can be adjusted. The larger
		// the better (probably). The motivation for this is to alleviate 
		// memory pressure, since you could theoretically want to follow
		// millions or billions of particles, and allocating them all at once
		// is wasteful and could cause you to run out of memory.
		constexpr int slot_cap = 131072;
		Slots::Slots slots {slot_cap};

		// All particles assumed to have same mass.
		slots.set_mass(opts.imp_mass_amu() * Constants::amu_to_kg);

		// Debug while setting things up, set all to dead to test filling
		std::cout << "setting all to dead\n";
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
#else

		// Turn off GPU option if it was turned on for a CPU-only build.
		if (opts.use_gpu_int() > 0)
		{
			std::cout << "GPU acceleration requested (use_gpu = " 
				<< opts.use_gpu() << ") but Flan was built without GPU "
				"support. Turning off.";
			opts.set_use_gpu("off");
		}

#endif
		
#ifdef USE_CUDA

		if (opts.use_gpu_int() > 0)
		{

			// Create device-side structs on each GPU if using GPUs
			std::vector<Slots::SlotsDevice> gpu_slots;
			std::vector<ImpurityStats::StatisticsDevice> gpu_stats;
			std::vector<Background::BackgroundDevice> gpu_bkgs;
			for (int dev = 0; dev < num_gpus; dev++) 
			{
				gpu_slots.push_back(slots.to_device(dev));
				gpu_stats.push_back(imp_stats.to_device(dev));
				gpu_bkgs.push_back(bkg.to_device(dev));
			}

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
					std::ref(gpu_stats[dev]), std::ref(oa_ioniz), 
					std::ref(oa_recomb), std::ref(opts), std::ref(timer));
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

				// Free memory of device-side structs
				Slots::free_slots(gpu_slots[dev], dev);
				Impurity::free_stats(gpu_stats[dev], dev);
				Background::free_bkg(gpu_bkgs[dev], dev);
			}

			return imp_stats;
		}

#endif

		// Dummy SlotsDevice, StatsDevice and BackgroundDevice. They're not 
		// used in CPU-only simulations.
		Slots::SlotsDevice slots_d {};
		ImpurityStats::StatisticsDevice imp_stats_d {};
		Background::BackgroundDevice bkg_d {};

		main_loop(slots, slots_d, bkg, bkg_d, imp_stats, imp_stats_d, oa_ioniz, 
			oa_recomb, opts, timer);

		return imp_stats;
	}
}
