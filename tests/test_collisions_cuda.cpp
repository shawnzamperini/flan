/*
* The GPU Nanbu collisions (collision_gpu) against the CPU's (nanbu_coll),
* particle by particle, with the thermal force on:
*
* - the slab with ion and electron heat fluxes along B, so every
*   background particle is sampled from the skewed distribution;
* - the HERMES test grid, whose heat fluxes are calculated from the
*   conduction and are interpolated along the field lines, over six zones.
*
* Slot i's generator on the GPU is pcg32(seed + i, i) (init_rngs_cuda), so
* the CPU is run one particle at a time with that generator. The two then
* draw the same random numbers, and the velocities must agree to rounding
* (the GPU fuses multiply-adds), sub-steps and rejected draws included. The
* skewed sampler's counts from the device must equal the CPU's exactly.
*/

#include <array>
#include <cmath>
#include <cuda_runtime.h>
#include <filesystem>
#include <gtest/gtest.h>
#include <omp.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "background.h"
#include "background_device.h"
#include "collisions.h"
#include "collisions_hd.h"
#include "impurity_stats.h"
#include "impurity_stats_device.h"
#include "multizone_cases.h"
#include "options.h"
#include "options_device.h"
#include "pcg32.h"
#include "read_bkg.h"
#include "slots.h"
#include "slots_device.h"

#include "collisions.cuh"
#include "init_rngs.cuh"

namespace
{
	constexpr std::uint64_t seed {91};

	/// Velocities after one collision step, and the sampler's counts
	struct Result
	{
		std::vector<double> vX;
		std::vector<double> vY;
		std::vector<double> vZ;
		Collisions::SkewCounts counts;
	};

	/// nanbu_coll on one particle at a time, each with slot i's GPU
	/// generator (every OpenMP thread gets a copy, since any may run it)
	Result cpu_one_by_one(const Slots::Slots& in,
		const Background::Background& bkg, const Options::Options& opts,
		const bool elec, const double dt)
	{
		Result r {};
		Collisions::reset_skew_counts();
		Impurity::Statistics stats {};
		for (int i {}; i < in.N(); ++i)
		{
			Slots::Slots one(1);
			one.set_mass(in.mass());
			one.set_Z(in.Z());
			one.set_state(0, in.state()[i]);
			one.set_zone(0, in.zone()[i]);
			one.set_q(0, in.q()[i]);
			one.set_t(0, in.t()[i]);
			one.set_x(0, in.x()[i]);
			one.set_y(0, in.y()[i]);
			one.set_z(0, in.z()[i]);
			one.set_vX(0, in.vX()[i]);
			one.set_vY(0, in.vY()[i]);
			one.set_vZ(0, in.vZ()[i]);
			std::vector<pcg32> rngs(omp_get_max_threads(),
				pcg32(seed + i, i));
			Collisions::nanbu_coll(one, bkg, opts, elec, dt, stats, rngs);
			r.vX.push_back(one.vX()[0]);
			r.vY.push_back(one.vY()[0]);
			r.vZ.push_back(one.vZ()[0]);
		}
		r.counts = Collisions::get_skew_counts(elec);
		Collisions::reset_skew_counts();
		return r;
	}

	/// collision_gpu on all of them, with the sampler's counts
	Result gpu_all(const Slots::Slots& in, Background::Background& bkg,
		Options::Options& opts, const bool elec, const double dt)
	{
		Background::BackgroundDevice bkg_d {bkg.to_device()};
		Slots::Slots gpu {in};
		Slots::SlotsDevice slots_d {gpu.to_device()};
		pcg32* rngs_d {init_rngs_cuda(slots_d, seed)};
		Options::OptionsDevice* opts_d {opts.to_device()};
		ImpurityStats::StatisticsDevice stats_d {};
		unsigned long long* skew_d {};
		const std::size_t skew_size {2 * Collisions::num_skew_counts
			* sizeof(unsigned long long)};
		cudaMalloc(&skew_d, skew_size);
		cudaMemset(skew_d, 0, skew_size);

		Collisions::collision_gpu(slots_d, bkg_d, elec, dt, stats_d, rngs_d,
			bkg.get_elec_mass_amu(), bkg.get_ion_mass_amu(), opts_d, skew_d);
		gpu = gpu.to_host(slots_d);
		const auto counts {Collisions::skew_counts_from_device(skew_d, 0)};

		cudaFree(skew_d);
		Options::free_opts(opts_d);
		free_rngs(rngs_d, slots_d.device_id);
		Slots::free_slots(slots_d);
		Background::free_bkg(bkg_d);

		Result r {};
		r.vX = gpu.vX();
		r.vY = gpu.vY();
		r.vZ = gpu.vZ();
		r.counts = counts[elec ? 1 : 0];
		return r;
	}

	/**
	* Runs both, for ions and electrons, and compares them.
	*
	* @param bkg Background (with heat fluxes).
	* @param opts Options (thermal force on).
	* @param slots The particles.
	* @param dt Time step (s).
	*/
	void collisions_match_cpu(Background::Background& bkg,
		Options::Options& opts, const Slots::Slots& slots, const double dt)
	{
		for (const bool elec : {false, true})
		{
			SCOPED_TRACE(elec ? "electrons" : "ions");
			const Result cpu {cpu_one_by_one(slots, bkg, opts, elec, dt)};
			const Result gpu {gpu_all(slots, bkg, opts, elec, dt)};

			// Relative to the particle's speed after the collision (or
			// before it, if that's larger), since a component can pass
			// near 0
			int bad {};
			for (int i {}; i < slots.N(); ++i)
			{
				const double v_in {std::hypot(slots.vX()[i], slots.vY()[i],
					slots.vZ()[i])};
				const double v_out {std::hypot(cpu.vX[i], cpu.vY[i],
					cpu.vZ[i])};
				const double scale {std::max(v_in, v_out)};
				const double diff {std::max({std::abs(gpu.vX[i] - cpu.vX[i]),
					std::abs(gpu.vY[i] - cpu.vY[i]),
					std::abs(gpu.vZ[i] - cpu.vZ[i])})};
				if (!(diff <= 1e-9 * scale) && bad++ < 5)
				{
					ADD_FAILURE() << "particle " << i << " (zone "
						<< slots.zone()[i] << ", q " << slots.q()[i]
						<< "): GPU (" << gpu.vX[i] << ", " << gpu.vY[i]
						<< ", " << gpu.vZ[i] << "), CPU (" << cpu.vX[i]
						<< ", " << cpu.vY[i] << ", " << cpu.vZ[i] << ")";
				}
			}
			EXPECT_EQ(bad, 0);

			// The skew was sampled, and the same way
			EXPECT_GT(cpu.counts.samples, 0);
			EXPECT_GT(cpu.counts.draws, cpu.counts.samples);
			EXPECT_EQ(gpu.counts.samples, cpu.counts.samples);
			EXPECT_EQ(gpu.counts.draws, cpu.counts.draws);
			EXPECT_EQ(gpu.counts.negative, cpu.counts.negative);
			EXPECT_EQ(gpu.counts.eps_capped, cpu.counts.eps_capped);
		}
	}
}

// Uniform slab with strong heat fluxes along B (eps 0.1 and above), and a
// time step long enough that the electrons are sub-cycled
TEST(CollisionsCUDA, MatchesCPUOnSlabWithHeatFlux)
{
	Options::Options opts {};
	opts.set_bkg_source("slab");
	opts.set_slab_ne(1e19);
	opts.set_slab_te(10.0);
	opts.set_slab_ti(10.0);
	opts.set_slab_bmag(1.0);
	opts.set_slab_heat_flux_i(0.3);
	opts.set_slab_heat_flux_e(0.3);
	Background::Background bkg {Background::read_bkg(opts)};
	ASSERT_FALSE(bkg.get_q_par_i().get_data().empty());

	Slots::Slots slots {MultiZone::spread_particles(bkg, 4000, 3)};
	collisions_match_cpu(bkg, opts, slots, 1e-7);
}

// HERMES test grid: heat fluxes from the conduction, interpolated along
// the field lines, in every zone
TEST(CollisionsCUDA, MatchesCPUOnHermesGrid)
{
	const auto dir {std::filesystem::temp_directory_path()
		/ ("flan_test_collisions_cuda_" + std::to_string(getpid()))};
	Background::Background bkg {MultiZone::hermes_background(dir)};
	std::filesystem::remove_all(dir);
	ASSERT_FALSE(bkg.get_q_par_e().get_data().empty());

	Options::Options opts {};
	Slots::Slots slots {MultiZone::spread_particles(bkg, 500, 5)};
	collisions_match_cpu(bkg, opts, slots, 1e-8);
}
