#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <thread>
#include <tuple>
#include <vector>

#ifdef USE_CUDA
#include <cuda_runtime.h>
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
#include "random.h"
#include "timer.h"
#include "utilities.h"
#include "variance_reduction.h"
#include "vectors.h"



namespace ImpurityTransport
{
	// Slots that hold a fixed number of particles and their properties
	struct ParticleSlots 
	{
		int N;        // total number of slots
		double* x;
		double* y;
		double* z;
		double* vx;
		double* vy;
		double* vz;
		double* vX;
		double* vY;
		double* vZ;
		int* alive;    // 1 = alive, 0 = dead
		int use_gpu;   // 0 = no, 1 = CUDA
		int device_id;
	};

	// Create empty slots to hold a batch of particles
	ParticleSlots create_slots(int N, const Options::Options& opts, 
		int device_id = 0) 
	{
		ParticleSlots slots;
		slots.N = N;
		slots.use_gpu = opts.use_gpu_int();
		slots.device_id = device_id;

		// CPU allocation
		if (opts.use_gpu_int() == 0) 
		{
			// Zero-initialize all arrays. If you add new arrays, don't forget 
			// to free them in free_slots!
			slots.x = new double[N]();
			slots.y = new double[N]();
			slots.z = new double[N]();
			slots.vx = new double[N]();
			slots.vy = new double[N]();
			slots.vz = new double[N]();
			slots.vX = new double[N]();
			slots.vY = new double[N]();
			slots.vZ = new double[N]();
			slots.alive = new int[N]();  // 0, all start dead
		}

#ifdef USE_CUDA

		// GPU allocation
		cudaSetDevice(device_id);

		cudaMalloc(&slots.x, N * sizeof(double));
		cudaMalloc(&slots.y, N * sizeof(double));
		cudaMalloc(&slots.z, N * sizeof(double));
		cudaMalloc(&slots.vx, N * sizeof(double));
		cudaMalloc(&slots.vy, N * sizeof(double));
		cudaMalloc(&slots.vz, N * sizeof(double));
		cudaMalloc(&slots.vX, N * sizeof(double));
		cudaMalloc(&slots.vY, N * sizeof(double));
		cudaMalloc(&slots.vZ, N * sizeof(double));
		cudaMalloc(&slots.alive, N * sizeof(int));

		// Initialize to zero
		cudaMemset(slots.x, 0, N * sizeof(double));
		cudaMemset(slots.y, 0, N * sizeof(double));
		cudaMemset(slots.z, 0, N * sizeof(double));
		cudaMemset(slots.vx, 0, N * sizeof(double));
		cudaMemset(slots.vy, 0, N * sizeof(double));
		cudaMemset(slots.vz, 0, N * sizeof(double));
		cudaMemset(slots.vX, 0, N * sizeof(double));
		cudaMemset(slots.vY, 0, N * sizeof(double));
		cudaMemset(slots.vZ, 0, N * sizeof(double));
		cudaMemset(slots.alive, 0, N * sizeof(int));

#endif
		return slots;

	}
	
	// Free ParticleSlots
	void free_slots(ParticleSlots& slots, int device_id = 0) 
	{
		if (!slots.use_gpu) 
		{
			delete[] slots.x;
			delete[] slots.y;
			delete[] slots.z;
			delete[] slots.vx;
			delete[] slots.vy;
			delete[] slots.vz;
			delete[] slots.vX;
			delete[] slots.vY;
			delete[] slots.vZ;
			delete[] slots.alive;
			return;
		}

#ifdef USE_CUDA

		cudaSetDevice(device_id);

		cudaFree(slots.x);
		cudaFree(slots.y);
		cudaFree(slots.z);
		cudaFree(slots.vx);
		cudaFree(slots.vy);
		cudaFree(slots.vz);
		cudaFree(slots.vX);
		cudaFree(slots.vY);
		cudaFree(slots.vZ);
		cudaFree(slots.alive);

#endif

	}

	void fill_slots(ParticleSlots& slots, int& rem_parts)
	{
	}

	
	void main_loop(ParticleSlots& slots, const Background::Background& bkg, 
		Impurity::Statistics& imp_stats, const OpenADAS::OpenADAS& oa_ioniz, 
		const OpenADAS::OpenADAS& oa_recomb, Options::Options& opts, 
		Timer::Timer& timer)
	{
		// Tracker for number of remaining particles to follow
		int rem_parts {opts.imp_num()};
		
		// Initial fill of slots with particles. 
		fill_slots(slots, rem_parts);

	
	}

#ifdef USE_CUDA

	// Function to enable each thread to control each GPU independently to 
	// allow them to run concurrently.
	void gpu_worker(ParticleSlots& slots, const Background::Background& bkg,
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
		// the better (probably). The collection of all slots is called the
		// casino baby cha-ching.
		constexpr int slot_cap = 131072;
		std::vector<ParticleSlots> casino;

#ifdef USE_CUDA

		// Get number of GPUs
		int num_gpus {};
		cudaGetDeviceCount(&num_gpus);
		std::cout << "Using " << num_gpus << " GPU(s)\n";

		// Create slots to hold particles on each GPU
		for (int dev = 0; dev < num_gpus; dev++) 
		{
			casino.push_back(create_slots(slot_cap, opts, dev));
		}

		// Spawn threads to launch main_loop on each available GPU
		std::vector<std::thread> threads;
		for (int dev = 0; dev < num_gpus; dev++) 
		{
			// Start a GPU worker on each thread. std::ref is used here 
			// because std::thread copies by default, which isn't necessary
			// or desirable for these larger objects.
			threads.emplace_back(gpu_worker, std::ref(casino[dev]),
				std::ref(bkg), std::ref(imp_stats), std::ref(oa_ioniz),
				std::ref(oa_recomb), std::ref(opts), std::ref(timer));
		}

		// Wait for all threads to finish
		for (auto& t : threads) {
			t.join();
		}
	
#else
		// CPUs don't need all the mumbo jumbo that GPUs need to utilize
		// multiple of them. With CPUs it's handled by appropriately assigning
		// one MPI task per CPU and then setting the correct number of OpenMP
		// tasks per CPU. That's used in main_loop when running on CPUs.
		ParticleSlots slots {create_slots(slot_cap, opts)};
		main_loop(slots, bkg, imp_stats, oa_ioniz, oa_recomb, opts, timer);
#endif

		// Free memory
		for (auto& s : casino) {
			free_slots(s, s.device_id);
		}



	
		return imp_stats;
	}
}
