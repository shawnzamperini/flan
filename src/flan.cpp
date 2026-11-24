/**
* @file  flan.cpp
*
* @brief Entry point and top level control of Flan library.
*
* An input file builds an Inputs map and calls flan with it, between
* MPI_Init and MPI_Finalize (see regression/blank_input_file.cpp).
*/
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <optional>
#include <string>

#include "background.h"
#include "boundary.h"
#include "collisions.h"
#include "config.h"
#include "errors.h"
#include "flan.h"
#include "flan_types.h"
#include "impurity_stats.h"
#include "impurity_transport.h"
#include "load_options.h"
#include "mpi.h"
#include "options.h"
#include "read_bkg.h"
#include "read_gkyl.h"
#include "save_results.h"
#include "sheath.h"
#include "sputtering.h"
#include "timer.h"
#include "vectors.h"

/**
* Runs a Flan simulation with the given input options.
*
* The steps are:
*
* 1. Load the options from inpts, filling in defaults for any not given
*    (Options::load_options).
* 2. Read the background plasma on rank 0 (Background::read_bkg) and
*    broadcast it to every rank.
* 3. With imp_source_opt = "sputter", find the wall surfaces and the
*    rate the background sputters atoms off them (Sputtering::Source),
*    on every rank. The particles then start from it, each standing for
*    its total rate / imp_num atoms/s, and where they land back on the
*    wall is tallied (Sputtering::Deposition).
* 4. Follow the impurities on every rank
*    (ImpurityTransport::follow_impurities).
* 5. Combine the statistics from every rank onto rank 0
*    (Impurity::reduce_stats), and calculate the density, velocities,
*    charge and Nanbu s from them.
* 6. Save the results on rank 0 (SaveResults::save_results) and print the
*    timing summary.
*
* MPI must already be initialized when flan is called, and the caller
* finalizes it afterwards. Every rank must call flan.
*
* @param inpts Input options from the input file. Options not in inpts keep
*        their default values.
*/
void flan(const Inputs& inpts)
{
	// MPI rank and number of processes
	int rank {};
	int nprocs {};
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// Create Timer object and begin total simulation timer
	Timer::Timer timer {};
	timer.start(Timer::Section::Total);

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

	// Without the thermal force, Nanbu collides impurities with
	// Maxwellians, which have no heat flux and so give no thermal force.
	// (With it on but no heat fluxes, check_thermal_force says so.)
	if (rank == 0 && opts.imp_collisions_int() == 1
		&& opts.imp_thermal_force_int() != 1)
	{
		std::cout << "\nWarning! imp_thermal_force = \"off\", so the Nanbu "
			"collision model\ncollides impurities with Maxwellians of the "
			"background plasma,\nwhich leaves out the thermal forces (from "
			"the temperature\ngradients). Set imp_thermal_force = \"on\" "
			"with a SOLEDGE3X\nbackground whose files have the face fluxes "
			"to include them.\n\n";
	}

	// Load plasma background on rank 0. The other ranks get it from the
	// broadcast below.
	timer.start(Timer::Section::Read);
	Background::Background bkg {};
	if (rank == 0) bkg = Background::read_bkg(opts);
	timer.stop(Timer::Section::Read);

	// Interpolate additional frames between each Gkeyll frame to artificially
	// increase the time resolution of the simulation.
	// To-do: Implement this later once I have an easier way to test it.

	// Broadcast the background to all the processes.
	bkg.broadcast(MPI_COMM_WORLD);

	// The thermal force needs the heat fluxes and the CPU
	Collisions::check_thermal_force(bkg, opts);

	// The charge limit for tests is only applied on the CPU
	if (opts.imp_iz_max_charge() >= 0 && opts.use_gpu_int() > 0)
	{
		Errors::fatal_error("imp_iz_max_charge only works on the CPU. Set "
			"use_gpu = \"off\".");
	}

	// The sputtered source. Every rank calculates its own from the
	// broadcast background, since every rank starts particles from it.
	// Each particle then stands for source_rate / imp_num atoms/s (its
	// weight), in place of imp_source_scale_fact / imp_num.
	std::optional<Sputtering::Source> source {};
	std::optional<Sputtering::Deposition> deposition {};
	if (opts.imp_source_opt_int() == 1)
	{
		if (opts.sput_table_file() == "undefined")
		{
			Errors::fatal_error("imp_source_opt = \"sputter\" needs a "
				"sputtering yield table, set with sput_table_file.");
		}
		if (opts.use_gpu_int() > 0)
		{
			Errors::fatal_error("imp_source_opt = \"sputter\" only runs on "
				"the CPU for now. Set use_gpu = \"off\".");
		}
		if (rank == 0) std::cout << "Calculating the sputtered source...\n";
		const Sputtering::Surface surface(bkg, opts.sput_max_wall_dist());
		const Sputtering::YieldTable table(opts.sput_table_file());

		// The ions' impact energies and angles from D+ test ions pushed
		// through the sheath field, which must then be on, or from the
		// formula at sput_angle_val
		Sheath::Model impact_model {Sheath::Model::off};
		if (opts.sput_impact_opt_int() == 1)
		{
			if (opts.sheath_model_int() == 0)
			{
				Errors::fatal_error("sput_impact_opt = \"sheath\" pushes the "
					"ions through the sheath field, but sheath_model is "
					"\"off\". Turn the sheath on, or use sput_impact_opt = "
					"\"formula\".");
			}
			impact_model = static_cast<Sheath::Model>(
				opts.sheath_model_int());
			if (rank == 0)
			{
				std::cout << "Pushing D+ test ions through the sheath for "
					"their impact energies and angles...\n";
			}
		}
		source.emplace(surface, bkg, table, opts, impact_model);
		if (rank == 0)
		{
			surface.print_summary(std::cout);
			source->print_summary(std::cout);
		}
		if (!(source->source_rate() > 0.0))
		{
			Errors::fatal_error("Nothing is sputtered off the wall (the "
				"source is 0 everywhere), so there are no particles to "
				"follow.");
		}
		// Where the sputtered atoms land back on the wall
		deposition.emplace(source->surface(), opts.imp_atom_num(),
			bkg.get_phi_sign());
		if (rank == 0 && opts.imp_source_scale_fact() != 1.0)
		{
			std::cout << "Warning! imp_source_scale_fact isn't used with "
				"imp_source_opt = \"sputter\".\nThe sputtered source ("
				<< source->source_rate() << " atoms/s) is used instead.\n";
		}
	}

	// The sheath field in front of the walls, from the same wall surface
	// as the sputtered source (or its own if there's no source). It's on
	// by default, so a background without walls just goes without, and on
	// the GPU (where it isn't written yet) it's left out with a warning.
	std::optional<Sheath::Sheath> sheath {};
	if (opts.sheath_model_int() > 0 && (source
		|| Sputtering::has_surfaces(bkg)))
	{
		if (opts.use_gpu_int() > 0)
		{
			if (rank == 0)
			{
				std::cout << "Warning! The sheath field (sheath_model = \""
					<< opts.sheath_model() << "\") only runs on the CPU for "
					"now, so this GPU run goes without it.\n";
			}
		}
		else
		{
			const auto model {static_cast<Sheath::Model>(
				opts.sheath_model_int())};
			if (source) sheath.emplace(source->surface(), bkg, model,
				opts.sheath_fd(), opts.sheath_drop());
			else
			{
				const Sputtering::Surface surface(bkg,
					opts.sput_max_wall_dist());
				sheath.emplace(surface, bkg, model, opts.sheath_fd(),
					opts.sheath_drop());
			}
			if (rank == 0) sheath->print_summary(std::cout);
		}
	}

	// Self-sputtering: impurities landing on the wall sputter more atoms
	// and carry on as one of them. With the sheath field on, they've been
	// accelerated through it already; without, the drop is added to their
	// impact energy.
	std::optional<Sputtering::SelfSource> self {};
	if (opts.sput_self_int() == 1)
	{
		if (!source)
		{
			Errors::fatal_error("sput_self = \"on\" needs the sputtered "
				"source (imp_source_opt = \"sputter\").");
		}
		if (opts.sput_self_table_file() == "undefined")
		{
			Errors::fatal_error("sput_self = \"on\" needs a "
				"self-sputtering yield table, set with "
				"sput_self_table_file.");
		}
		const Sputtering::YieldTable self_table(opts.sput_self_table_file());
		self.emplace(*source, self_table, opts, sheath.has_value());
		deposition->set_self(&*self);
		if (rank == 0)
		{
			std::cout << "Self-sputtering on, yields from "
				<< opts.sput_self_table_file() << " (" << self_table
				.energy_min() << " to " << self_table.energy_max()
				<< " eV). Roulette below " << self->min_weight()
				<< " atoms/s; the run stops at generation "
				<< opts.sput_self_max_gen() << ".\n";
		}
	}

	// Begin main particle following loop.
	timer.start(Timer::Section::Imp);
	Boundary::Fates fates {};
	Collisions::reset_skew_counts();
	Impurity::Statistics imp_stats {ImpurityTransport::follow_impurities(bkg, opts, 
		timer, source ? &*source : nullptr,
		deposition ? &*deposition : nullptr, &fates,
		sheath ? &*sheath : nullptr)};
	timer.stop(Timer::Section::Imp);

	// How the thermal force's skewed sampling went, over every rank
	std::optional<std::array<Collisions::SkewCounts, 2>> skew_counts {};
	if (opts.imp_thermal_force_int() == 1)
	{
		skew_counts = Collisions::reduce_skew_counts();
		if (rank == 0) Collisions::print_skew_counts(*skew_counts);
	}

	// Add up where the sputtered atoms landed over every rank
	if (deposition)
	{
		deposition->reduce(0);
		if (rank == 0)
			deposition->print_summary(std::cout, source->source_rate());
	}
	if (self)
	{
		self->reduce(0);
		if (rank == 0) self->print_summary(std::cout, source->source_rate());
	}

	// Reduce imp_stats from each rank, defined in impurity_stats.cpp. Only
	// rank 0 returns a meaningful result here.
	Impurity::Statistics global_stats {Impurity::reduce_stats(imp_stats)};

	// Convert the statistics into meaningful quantities.
	if (rank == 0)
	{	
		std::cout << "Calculating derived quantities...\n";
		std::cout << "  Density...\n";
		global_stats.calc_density(bkg, opts);
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
		timer.start(Timer::Section::Save);
		SaveResults::save_results(bkg, global_stats, opts,
			source ? &*source : nullptr,
			deposition ? &*deposition : nullptr, &fates,
			skew_counts ? &*skew_counts : nullptr);
		timer.stop(Timer::Section::Save);
	}

	// End total timer, then print summary just for MPI rank 0. It should be
	// about the same for all the other ranks. This is also going to give
	// you the wall time effectively, which is what we're after.
	timer.stop(Timer::Section::Total);
	if (rank == 0)
	{
		timer.print_summary();
	}

}  // flan
