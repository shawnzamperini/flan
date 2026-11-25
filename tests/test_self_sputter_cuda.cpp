/*
* Self-sputtering on the GPU against the CPU, on the multi-zone test
* backgrounds:
*
* - the HERMES test grid's targets (wall faces), without the sheath field,
*   so an ion's impact energy gets the floating drop;
* - the SOLEDGE3X test mesh's staircase (atoms absorbed in wall cells),
*   with the sheath field on (no drop added).
*
* Tungsten particles are put past the faces of their zones, as in the
* deposition tests (test_sputtering_cuda), with weights around the source
* particles' so that roulette is played, and generations below the
* limit. The boundary step (check_bounds_gpu with a SelfDevice) is
* compared with check_bounds_cpu run one particle at a time with slot i's
* GPU generator, pcg32(seed + i, i): every slot's state afterward (the new
* atom of those that carry on), the self-sputtering tallies, the
* deposition tallies and the count of self-sputtered particles. They draw
* the same random numbers, so the slots must agree to rounding and the
* tallies, added in a different order, to 1e-12.
*
* A last test puts every particle at the generation limit and checks
* that the GPU raises the avalanche flag, which stops the run, only then.
*/

#include <algorithm>
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
#include "boundary.h"
#include "boundary_hd.h"
#include "multizone_cases.h"
#include "options.h"
#include "pcg32.h"
#include "slots.h"
#include "slots_device.h"
#include "sputtering.h"
#include "sputtering_hd.h"

#include "boundary.cuh"
#include "init_rngs.cuh"

namespace
{
	constexpr std::uint64_t seed {31};
	constexpr int max_gen {6};
	constexpr int imp_num {1000};

	const std::string d_on_w {std::string {FLAN_SOURCE_DIR}
		+ "/data/sputtering/D_on_W.h5"};
	const std::string w_on_w {std::string {FLAN_SOURCE_DIR}
		+ "/data/sputtering/W_on_W.h5"};

	std::filesystem::path temp_dir(const std::string& name)
	{
		return std::filesystem::temp_directory_path() / ("flan_test_"
			+ name + "_" + std::to_string(getpid()));
	}

	Options::Options self_options()
	{
		Options::Options opts {};
		opts.set_imp_atom_num(74);
		opts.set_imp_mass_amu(183.84);
		opts.set_imp_num(imp_num);
		opts.set_sput_self_min_weight(0.5);
		opts.set_sput_self_max_gen(max_gen);
		return opts;
	}

	double rel(const double a, const double b, const double scale)
	{
		return std::abs(a - b) / std::max({std::abs(a), std::abs(b),
			scale});
	}

	/// Particles past the faces of their zones, with weights from a tenth
	/// to twice the source particles' w0, and generation below the limit
	/// (or at it, if at_limit)
	Slots::Slots landing_particles(const Background::Background& bkg,
		const int per_zone, const double w0, const bool at_limit)
	{
		Slots::Slots slots {MultiZone::spread_particles(bkg, per_zone, 23)};
		pcg32 rng(9, 4);
		const auto& zones {bkg.get_zones()};
		for (int i {}; i < slots.N(); ++i)
		{
			const int zone {slots.zone()[i]};
			const double pos[3] {slots.x()[i], slots.y()[i], slots.z()[i]};
			double moved[3] {};
			for (int d {}; d < 3; ++d)
			{
				const double span {zones.get_max(zone, d)
					- zones.get_min(zone, d)};
				moved[d] = pos[d] + 0.6 * span * (rng.next_double() - 0.5);
			}
			slots.set_x(i, moved[0]);
			slots.set_y(i, moved[1]);
			slots.set_z(i, moved[2]);
			slots.set_pt(i, std::pow(10.0, -9.0 + 5.0 * rng.next_double()));
			slots.set_prompt_by(i, rng.next_double() < 0.5 ? -1.0
				: slots.pt()[i] * (0.5 + rng.next_double()));
			slots.set_weight(i, w0 * (0.1 + 1.9 * rng.next_double()));
			slots.set_generation(i, at_limit ? max_gen
				: static_cast<int>(max_gen * rng.next_double()) % max_gen);
			// From a fifth to 4 times spread_particles' 1e4 m/s, so some
			// neutrals land below the W on W threshold and some above
			const double f {0.2 + 3.8 * rng.next_double()};
			slots.set_vX(i, f * slots.vX()[i]);
			slots.set_vY(i, f * slots.vY()[i]);
			slots.set_vZ(i, f * slots.vZ()[i]);
			if (i % 5 == 0) slots.set_q(i, 0);
		}
		return slots;
	}

	/// A copy of slot i of slots, alone
	Slots::Slots one_slot(const Slots::Slots& slots, const int i)
	{
		Slots::Slots one(1);
		one.set_mass(slots.mass());
		one.set_Z(slots.Z());
		one.set_state(0, slots.state()[i]);
		one.set_zone(0, slots.zone()[i]);
		one.set_t(0, slots.t()[i]);
		one.set_tidx(0, slots.tidx()[i]);
		one.set_x(0, slots.x()[i]);
		one.set_y(0, slots.y()[i]);
		one.set_z(0, slots.z()[i]);
		one.set_xidx(0, slots.xidx()[i]);
		one.set_yidx(0, slots.yidx()[i]);
		one.set_zidx(0, slots.zidx()[i]);
		one.set_vX(0, slots.vX()[i]);
		one.set_vY(0, slots.vY()[i]);
		one.set_vZ(0, slots.vZ()[i]);
		one.set_vx(0, slots.vx()[i]);
		one.set_vy(0, slots.vy()[i]);
		one.set_vz(0, slots.vz()[i]);
		one.set_q(0, slots.q()[i]);
		one.set_pt(0, slots.pt()[i]);
		one.set_prompt_by(0, slots.prompt_by()[i]);
		one.set_weight(0, slots.weight()[i]);
		one.set_generation(0, slots.generation()[i]);
		return one;
	}

	/**
	* Runs the boundary step with self-sputtering on both and compares
	* them.
	*
	* @return How many particles carried on as self-sputtered atoms (CPU).
	*/
	long long self_sputter_matches_cpu(Background::Background& bkg,
		const bool sheath_on)
	{
		const Options::Options opts {self_options()};
		const Sputtering::Surface surface(bkg);
		const Sputtering::YieldTable table(d_on_w);
		const Sputtering::YieldTable self_table(w_on_w);
		const Sputtering::Source source(surface, bkg, table, opts);
		const double w0 {source.source_rate() / imp_num};
		const Slots::Slots slots {landing_particles(bkg, 2000, w0, false)};

		// CPU, one at a time with slot i's generator
		Sputtering::Deposition cpu_dep(surface, 74, bkg.get_phi_sign());
		Sputtering::SelfSource cpu_self(source, self_table, opts, sheath_on);
		cpu_dep.set_self(&cpu_self);
		Boundary::Fates fates {};
		std::vector<Slots::Slots> cpu;
		cpu.reserve(slots.N());
		for (int i {}; i < slots.N(); ++i)
		{
			Slots::Slots one {one_slot(slots, i)};
			std::vector<pcg32> rngs(omp_get_max_threads(),
				pcg32(seed + i, i));
			Boundary::check_bounds_cpu(one, bkg, opts, rngs, fates, &cpu_dep);
			cpu.push_back(one);
		}

		// GPU
		Sputtering::Deposition gpu_dep(surface, 74, bkg.get_phi_sign());
		Sputtering::SelfSource gpu_self(source, self_table, opts, sheath_on);
		Background::BackgroundDevice bkg_d {bkg.to_device()};
		Slots::Slots copy {slots};
		Slots::SlotsDevice slots_d {copy.to_device()};
		pcg32* rngs_d {init_rngs_cuda(slots_d, seed)};
		Sputtering::DepositionDevice dep_d {gpu_dep.to_device(0)};
		Sputtering::SelfDevice self_d {gpu_self.to_device(0)};
		unsigned long long* fates_d {};
		cudaMalloc(&fates_d, Boundary::num_fate_indices
			* sizeof(unsigned long long));
		cudaMemset(fates_d, 0, Boundary::num_fate_indices
			* sizeof(unsigned long long));
		Boundary::check_bounds_gpu(slots_d, bkg_d,
			Boundary::make_params(bkg, opts, bkg_d.chi), rngs_d, fates_d,
			&dep_d, &self_d);
		EXPECT_EQ(cudaDeviceSynchronize(), cudaSuccess);
		EXPECT_FALSE(Sputtering::avalanche_on_device(self_d));
		const Slots::Slots gpu {copy.to_host(slots_d)};
		gpu_dep.add_from_device(dep_d);
		gpu_self.add_from_device(self_d);
		Boundary::Fates gpu_fates {};
		Boundary::add_fates_from_device(fates_d, 0, gpu_fates);
		cudaFree(fates_d);
		Sputtering::free_self(self_d);
		Sputtering::free_deposition(dep_d);
		free_rngs(rngs_d, slots_d.device_id);
		Slots::free_slots(slots_d);
		Background::free_bkg(bkg_d);

		// Each slot
		int bad {};
		for (int i {}; i < slots.N(); ++i)
		{
			const Slots::Slots& c {cpu[i]};
			const double v {std::max(1.0, std::hypot(c.vX()[0], c.vY()[0],
				c.vZ()[0]))};
			bool ok {gpu.state()[i] == c.state()[0]};
			ok &= gpu.zone()[i] == c.zone()[0];
			ok &= gpu.xidx()[i] == c.xidx()[0];
			ok &= gpu.yidx()[i] == c.yidx()[0];
			ok &= gpu.zidx()[i] == c.zidx()[0];
			ok &= gpu.q()[i] == c.q()[0];
			ok &= gpu.generation()[i] == c.generation()[0];
			ok &= rel(gpu.x()[i], c.x()[0], 1e-6) < 1e-12;
			ok &= rel(gpu.y()[i], c.y()[0], 1e-6) < 1e-12;
			ok &= rel(gpu.z()[i], c.z()[0], 1e-6) < 1e-12;
			ok &= rel(gpu.vX()[i], c.vX()[0], v) < 1e-10;
			ok &= rel(gpu.vY()[i], c.vY()[0], v) < 1e-10;
			ok &= rel(gpu.vZ()[i], c.vZ()[0], v) < 1e-10;
			ok &= rel(gpu.weight()[i], c.weight()[0], 1e-300) < 1e-12;
			ok &= gpu.pt()[i] == c.pt()[0];
			ok &= gpu.prompt_by()[i] == c.prompt_by()[0];
			ok &= gpu.t()[i] == c.t()[0];
			if (!ok && bad++ < 5)
			{
				ADD_FAILURE() << "slot " << i << ": GPU state "
					<< gpu.state()[i] << " gen " << gpu.generation()[i]
					<< " w " << gpu.weight()[i] << " x (" << gpu.x()[i]
					<< ", " << gpu.y()[i] << ", " << gpu.z()[i] << ") v ("
					<< gpu.vX()[i] << ", " << gpu.vY()[i] << ", "
					<< gpu.vZ()[i] << "); CPU state " << c.state()[0]
					<< " gen " << c.generation()[0] << " w "
					<< c.weight()[0] << " x (" << c.x()[0] << ", "
					<< c.y()[0] << ", " << c.z()[0] << ") v (" << c.vX()[0]
					<< ", " << c.vY()[0] << ", " << c.vZ()[0] << ")";
			}
		}
		EXPECT_EQ(bad, 0);

		// The tallies
		auto same = [](const std::vector<double>& g,
			const std::vector<double>& c, const char* name)
		{
			ASSERT_EQ(g.size(), c.size()) << name;
			int n_bad {};
			for (std::size_t n {}; n < c.size(); ++n)
			{
				if (!(rel(g[n], c[n], 1e-300) < 1e-12) && n_bad++ < 3)
				{
					ADD_FAILURE() << name << "[" << n << "]: GPU " << g[n]
						<< ", CPU " << c[n];
				}
			}
		};
		same(gpu_self.self_sputtered(), cpu_self.self_sputtered(),
			"self_sputtered");
		same(gpu_self.generation_weight(), cpu_self.generation_weight(),
			"generation_weight");
		same(gpu_self.generation_landed(), cpu_self.generation_landed(),
			"generation_landed");
		same(gpu_dep.redeposited(), cpu_dep.redeposited(), "redeposited");
		same(gpu_dep.charge(), cpu_dep.charge(), "charge");
		EXPECT_LT(rel(gpu_self.self_total(), cpu_self.self_total(), 1e-300),
			1e-12);
		EXPECT_LT(rel(gpu_self.landed(), cpu_self.landed(), 1e-300), 1e-12);
		EXPECT_LT(rel(gpu_self.mean_energy(), cpu_self.mean_energy(),
			1e-300), 1e-12);
		EXPECT_LT(rel(gpu_self.mean_angle(), cpu_self.mean_angle(), 1e-300),
			1e-12);
		EXPECT_LT(rel(gpu_self.roulette_lost(), cpu_self.roulette_lost(),
			1e-300), 1e-12);
		EXPECT_EQ(gpu_self.num_below_threshold(),
			cpu_self.num_below_threshold());
		EXPECT_EQ(gpu_self.num_above_table(), cpu_self.num_above_table());
		EXPECT_EQ(gpu_fates.self_sputtered, fates.self_sputtered);
		EXPECT_EQ(gpu_fates.wall, fates.wall);
		EXPECT_EQ(gpu_fates.chi, fates.chi);

		// Something of each kind happened
		EXPECT_GT(cpu_self.landed(), 0.0);
		EXPECT_GT(cpu_self.roulette_lost(), 0.0);
		EXPECT_GT(cpu_self.num_below_threshold(), 0);
		return fates.self_sputtered;
	}
}

// HERMES targets, the drop added to ions' impact energy
TEST(SelfSputterCUDA, MatchesCPUAtHermesTargets)
{
	if (!std::filesystem::exists(w_on_w)) GTEST_SKIP() << w_on_w;
	const auto dir {temp_dir("self_hermes")};
	Background::Background bkg {MultiZone::hermes_background(dir)};
	std::filesystem::remove_all(dir);
	EXPECT_GT(self_sputter_matches_cpu(bkg, false), 50);
}

// SOLEDGE3X staircase, absorbed in wall cells, sheath field on
TEST(SelfSputterCUDA, MatchesCPUInS3XWallCells)
{
	if (!std::filesystem::exists(w_on_w)) GTEST_SKIP() << w_on_w;
	const auto dir {temp_dir("self_s3x")};
	Background::Background bkg {MultiZone::s3x_background(dir)};
	std::filesystem::remove_all(dir);
	EXPECT_GT(self_sputter_matches_cpu(bkg, true), 50);
}

// Particles at sput_self_max_gen that would carry on raise the avalanche
// flag; the same particles a generation lower don't
TEST(SelfSputterCUDA, AvalancheFlag)
{
	if (!std::filesystem::exists(w_on_w)) GTEST_SKIP() << w_on_w;
	const auto dir {temp_dir("self_avalanche")};
	Background::Background bkg {MultiZone::hermes_background(dir)};
	std::filesystem::remove_all(dir);
	const Options::Options opts {self_options()};
	const Sputtering::Surface surface(bkg);
	const Sputtering::Source source(surface, bkg,
		Sputtering::YieldTable {d_on_w}, opts);
	const double w0 {source.source_rate() / imp_num};

	for (const bool at_limit : {false, true})
	{
		Sputtering::Deposition dep(surface, 74, bkg.get_phi_sign());
		Sputtering::SelfSource self(source, Sputtering::YieldTable {w_on_w},
			opts, false);
		Slots::Slots slots {landing_particles(bkg, 500, w0, at_limit)};
		if (!at_limit)
		{
			for (int i {}; i < slots.N(); ++i)
				slots.set_generation(i, max_gen - 1);
		}
		Background::BackgroundDevice bkg_d {bkg.to_device()};
		Slots::SlotsDevice slots_d {slots.to_device()};
		pcg32* rngs_d {init_rngs_cuda(slots_d, seed)};
		Sputtering::DepositionDevice dep_d {dep.to_device(0)};
		Sputtering::SelfDevice self_d {self.to_device(0)};
		Boundary::check_bounds_gpu(slots_d, bkg_d,
			Boundary::make_params(bkg, opts, bkg_d.chi), rngs_d, nullptr,
			&dep_d, &self_d);
		EXPECT_EQ(cudaDeviceSynchronize(), cudaSuccess);
		EXPECT_EQ(Sputtering::avalanche_on_device(self_d), at_limit);
		const auto st {Sputtering::self_stats_from_device(self_d)};
		EXPECT_GT(st.mean_energy, 0.0);
		Sputtering::free_self(self_d);
		Sputtering::free_deposition(dep_d);
		free_rngs(rngs_d, slots_d.device_id);
		Slots::free_slots(slots_d);
		Background::free_bkg(bkg_d);
	}
}
