#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <thread>
#include <tuple>
#include <vector>

#ifdef USE_CUDA
#include <cuda_runtime.h>
#include "particle_slots.cuh"
#endif

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
#include "particle_slots.h"
#include "random.h"
#include "timer.h"
#include "utilities.h"
#include "variance_reduction.h"
#include "vectors.h"



namespace ImpurityTransport
{



	
	void main_loop(ParticleSlots::Slots& slots, 
		const Background::Background& bkg, 
		Impurity::Statistics& imp_stats, const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADAS& oa_recomb, Options::Options& opts, 
		Timer::Timer& timer)
	{
		// Tracker for number of remaining particles to follow
		int rem_parts {opts.imp_num()};
		
		// Initial fill of slots with particles. 
		ParticleSlots::fill_slots(slots, rem_parts);

	
	}

#ifdef USE_CUDA

	// Function to enable each thread to control each GPU independently to 
	// allow them to run concurrently.
	void gpu_worker(ParticleSlots::Slots& slots, 
		const Background::Background& bkg,
        Impurity::Statistics& imp_stats, const OpenADAS::OpenADAS& oa_ioniz,
        const OpenADAS::OpenADAS& oa_recomb, Options::Options& opts,
        Timer::Timer& timer)
	{
		cudaSetDevice(slots.device_id);   // THIS THREAD USES THIS GPU
		main_loop(slots, bkg, imp_stats, oa_ioniz, oa_recomb, opts, timer);
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
		// Statistics object. Option to control if the three velocity arrays
		// are allocated (to save memory).
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

#ifdef USE_CUDA

		// Get number of GPUs
		int num_gpus {};
		cudaGetDeviceCount(&num_gpus);
		std::cout << "Using " << num_gpus << " GPU(s)\n";

		// Create slots to hold particles on each GPU
		std::vector<ParticleSlots::Slots> gpu_slots;
		for (int dev = 0; dev < num_gpus; dev++) 
		{
			gpu_slots.push_back(ParticleSlots::create_slots(slot_cap, opts.use_gpu_int(), dev));
		}

		// Spawn threads to launch main_loop on each available GPU
		std::vector<std::thread> threads;
		for (int dev = 0; dev < num_gpus; dev++) 
		{
			// Start a GPU worker on each thread. std::ref is used here 
			// because std::thread copies by default, which isn't necessary
			// or desirable for these larger objects.
			threads.emplace_back(gpu_worker, std::ref(gpu_slots[dev]),
				std::ref(bkg), std::ref(imp_stats), std::ref(oa_ioniz),
				std::ref(oa_recomb), std::ref(opts), std::ref(timer));
		}

		// Wait for all threads to finish
		for (auto& t : threads) {
			t.join();
		}

		// Free memory
		for (auto& s : gpu_slots) {
			ParticleSlots::free_slots(s, s.device_id);
		}
	
#else
		// CPUs don't need all the mumbo jumbo that GPUs need to utilize
		// multiple of them. With CPUs it's handled by appropriately assigning
		// one MPI task per CPU and then setting the correct number of OpenMP
		// tasks per CPU. That's used in main_loop when running on CPUs.
		ParticleSlots::Slots slots {ParticleSlots::create_slots(slot_cap, opts.use_gpu_int())};
		main_loop(slots, bkg, imp_stats, oa_ioniz, oa_recomb, opts, timer);

		// Free memory
		ParticleSlots::free_slots(slots);
#endif




	
		return imp_stats;
	}
}
