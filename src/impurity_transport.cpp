#include <algorithm>
#include <cmath>
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
#include "boundary.cuh"
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

	
	// Perform particle step for each particle in slots
	void step_cpu(Slots::Slots& slots, const Background::Background& bkg, 
		const Options::Options& opts, const double dt)
	{

		#pragma omp parallel for
		for (int i=0; i < slots.N(); ++i)
		{
			// Update time
			slots.set_t(i, slots.t()[i] + dt);

			// Update curvilinear position
			slots.set_x(i, slots.x()[i] + slots.vx()[i] * dt);
			slots.set_y(i, slots.y()[i] + slots.vy()[i] * dt);
			slots.set_z(i, slots.z()[i] + slots.vz()[i] * dt);
		}
	}

	// ------------------
	// Wrapper functions
	// ------------------

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
			// Defined in cuda/boris.cu
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
			// Defined in cuda/boundary.cu
			Boundary::check_bounds_gpu(slots_d, bkg_d, 
				opts.tbound_type_int(), 
				opts.min_xbound_type_int(), opts.max_xbound_type_int(),
				opts.min_ybound_type_int(), opts.max_ybound_type_int(),
				opts.min_zbound_type_int(), opts.max_zbound_type_int(),
				opts.imp_xbound_buffer(), opts.imp_ybound_buffer(),
				opts.imp_zbound_buffer(), opts.lcfs_x());
			return;
		}
#endif

		// Defined in boundary.cpp
		Boundary::check_bounds_cpu(slots, bkg, opts);
	}

	
	// Perform particle step
	void step_wrapper(Slots::Slots& slots, Slots::SlotsDevice& slots_d,
		const Background::Background& bkg, 
		const Background::BackgroundDevice& bkg_d, 
		const Options::Options& opts)
	{

#ifdef USE_CUDA
		if (opts.use_gpu_int() > 0) 
		{
			// Defined in impurity_transport.cu
			//step_gpu();
		}
#endif

		// Defined above
		step_cpu(slots, bkg, opts, opts.imp_time_step());

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

		std::cout << "initial\n";
		std::cout << " vx[0] = " << slots.vx()[0] << '\n';
		std::cout << "  t[0] = " << slots.t()[0] << '\n';
		std::cout << "  x[0] = " << slots.x()[0] << '\n';
		std::cout << "  y[0] = " << slots.y()[0] << '\n';
		std::cout << "  z[0] = " << slots.z()[0] << '\n';
		
		// Find starting grid index
		find_containing_cell_wrapper(slots, slots_d, bkg, bkg_d, opts);

		// Boris algorithm: The velocity stored in the Slots objects
		// will actually be the velocity at a half timestep earlier, i.e.,
		// at t - dt/2. So we still need to push the particle velocity
		// back by half a time step.
		boris_wrapper(slots, slots_d, bkg, bkg_d, opts, 
			-opts.imp_time_step() / 2.0);

		std::cout << "initial (after half step)\n";
		std::cout << " vx[0] = " << slots.vx()[0] << '\n';
		std::cout << "  t[0] = " << slots.t()[0] << '\n';
		std::cout << "  x[0] = " << slots.x()[0] << '\n';
		std::cout << "  y[0] = " << slots.y()[0] << '\n';
		std::cout << "  z[0] = " << slots.z()[0] << '\n';

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
			// To-do

			// We create a scope for each step so that we can use a scoped
			// timer to profile the time spent in each step of the loop.
			// A scoped timer is a bit safer because it automatically stops
			// when it goes out of scope.
			{
				Timer::ScopedTimer t(timer.acc(Timer::Section::Boris));

				// Update velocity (Boris).
				boris_wrapper(slots, slots_d, bkg, bkg_d, opts, 
					opts.imp_time_step());
			}

			// Perform particle step
			{
				Timer::ScopedTimer t(timer.acc(Timer::Section::Step));
				step_wrapper(slots, slots_d, bkg, bkg_d, opts);
			}

			// Bounds checking
			{
				Timer::ScopedTimer t(timer.acc(Timer::Section::Bounds));
				check_bounds_wrapper(slots, slots_d, bkg, bkg_d, opts);
			}

			// Update particle indices
			{
				Timer::ScopedTimer t(timer.acc(Timer::Section::FindCell));
				find_containing_cell_wrapper(slots, slots_d, bkg, bkg_d, opts);
			}

			// Record statistics
			{
				Timer::ScopedTimer t(timer.acc(Timer::Section::Record));
				record_stats_wrapper(imp_stats, imp_stats_d, slots, slots_d, 
					opts, opts.imp_time_step());
			}

			// Replace dead particles
			{
				Timer::ScopedTimer t(timer.acc(Timer::Section::FillSlots));
				fill_slots_wrapper(slots, slots_d, rem_parts, opts);
			}

			// User feedback
			// To-do

			// Check if all the particles in slots are dead. If this happens
			// after fill_slots, it means there were no more alive particles
			// to swap in and all the remaining ones are dead. So we're done.
			{
				Timer::ScopedTimer t(timer.acc(Timer::Section::AllDead));
				all_dead = all_dead_wrapper(slots, slots_d, opts, rem_parts);
			}

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

		// Slot capacity, 2^17 seems reasonable but can be adjusted. The larger
		// the better (probably). The motivation for this is to alleviate 
		// memory pressure, since you could theoretically want to follow
		// millions or billions of particles, and allocating them all at once
		// is wasteful and could cause you to run out of memory.
		constexpr int slot_cap = 131072;
		Slots::Slots slots {std::min(slot_cap, opts.imp_num())};

		// All particles assumed to have same mass.
		slots.set_mass(opts.imp_mass_amu() * Constants::amu_to_kg);

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

				// Initialize RNGs with SAME SEED FOR ALL OF THEM 
				init_slot_rngs(gpu_slots[dev], opts.seed());

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
