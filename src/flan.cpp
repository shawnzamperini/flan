/**
* @file  flan.cpp
*
* @brief Entry point and top level control of Flan library.
*/
#include <cmath>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <string>

#include "background.h"
#include "config.h"
#include "flan.h"
#include "flan_types.h"
#include "impurity_stats.h"
#include "impurity_transport.h"
#include "load_options.h"
#include "mpi.h"
#include "options.h"
#include "read_bkg.h"
#include "read_gkyl.h"
#include "read_test.h"
#include "save_results.h"
#include "timer.h"
#include "vectors.h"

void flan(const Inputs& inpts)
{
	// MPI rank and number of processes
	int rank {};
	int nprocs {};
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// Create Timer object and begin total simulation timer
	Timer::Timer timer {};
	timer.start_total_timer();

	// For printing output, flush buffer after every output operation
	std::cout << std::unitbuf;

	// Load version info from config.h file, print MPI/OMP usage
	if (rank == 0)
	{
		int omp_threads = 0;
		#pragma omp parallel
		{
			if (omp_get_thread_num() == 0)
				omp_threads = omp_get_num_threads();
		}

		std::cout <<
R"(
============================================================
               Flan (the code, not the dog)
============================================================
)" ;

		std::cout << " Version:           " << PROJECT_VERSION << "\n";
		std::cout << " Git Commit:        " << GIT_COMMIT_HASH << "\n";
		std::cout << " MPI Processes:     " << nprocs << "\n";
		std::cout << " OMP Threads/Rank:  " << omp_threads << "\n";

		std::cout <<
	"============================================================\n\n";
	}

	// Load options
	if (rank == 0) std::cout << "Loading input...\n";
	Options::Options opts {Options::load_options(inpts)};

	// Load plasma background. Returns a Background class option.
	timer.start_read_timer();
	Background::Background bkg {};
	if (rank == 0) bkg = Background::read_bkg(opts);
	timer.end_read_timer();

	// Interpolate additional frames between each Gkeyll frame to artificially
	// increase the time resolution of the simulation.
	// To-do: Implement this later once I have an easier way to test it.

	// Broadcast the background to all the processes.
	bkg.broadcast(MPI_COMM_WORLD);

	// Begin main particle following loop.
	timer.start_imp_timer();
	Impurity::Statistics imp_stats {Impurity::follow_impurities(bkg, opts, 
		timer)};
	timer.end_imp_timer();

	// Reduce imp_stats from each rank, defined in impurity_stats.cpp. Only
	// rank 0 returns a meaningful result here.
	Impurity::Statistics global_stats {Impurity::reduce_stats(imp_stats)};

	// Convert the statistics into meaningful quantities.
	if (rank == 0)
	{	
		std::cout << "Calculating derived quantities...\n";
		std::cout << "  Density...\n";
		global_stats.calc_density(bkg, opts.imp_num(), 
			opts.imp_source_scale_fact());
		std::cout << "  Velocity...\n";
		global_stats.calc_vels();
		std::cout << "  Charge...\n";
		global_stats.calc_charge();
		std::cout << "  Nanbu - s...\n";
		global_stats.calc_s();
	}

	// Only need to save data on root process
	if (rank == 0)
	{
		// Save simulation results.
		timer.start_save_timer();
		SaveResults::save_results(bkg, global_stats, opts);
		timer.end_save_timer();

		// End timer, print summary. Even though this is really just rank 0's
		// time summary, it should be about the same for every other rank 
		// (excluding the rank 0 specific tasks). It's still good for telling
		// you how long the program took in each section. It's essentially
		// the wall time.
		timer.end_total_timer();
		timer.print_summary();
	}

}
