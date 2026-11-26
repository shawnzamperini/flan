/*
* The sputtered source and the deposition on the GPU against the CPU, on
* the multi-zone test backgrounds:
*
* - births (fill_slots_gpu with a Sputtering::SourceDevice) against
*   Source::sample_birth: on the HERMES test grid (target faces) with the
*   impact energy from the formula, and on the SOLEDGE3X test mesh (a
*   staircase wall from chi) with the impact energies and angles of test
*   ions pushed through the sheath;
* - ionization with a charge limit (imp_iz_max_charge) and the prompt
*   window started on the first ionization, against OpenADAS::ioniz_recomb;
* - the deposition tallies from the boundary kernel against
*   check_bounds_cpu's, at the HERMES targets and in the SOLEDGE3X wall
*   cells.
*
* Slot i's generator on the GPU is pcg32(seed + i, i) (init_rngs_cuda), so
* the CPU is run one particle at a time with that generator. The two then
* draw the same random numbers and must agree to rounding. The tallies are
* sums added to in a different order on the GPU, so they agree to
* rounding too.
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
#include "multizone_cases.h"
#include "openadas.h"
#include "openadas_device.h"
#include "options.h"
#include "options_device.h"
#include "pcg32.h"
#include "read_bkg.h"
#include "sheath.h"
#include "slots.h"
#include "slots_device.h"
#include "sputtering.h"
#include "sputtering_hd.h"

#include "boundary.cuh"
#include "init_rngs.cuh"
#include "openadas.cuh"
#include "slots.cuh"

#include "openadas_test_utils.h"

namespace
{
	constexpr std::uint64_t seed {29};

	const std::string d_on_w {std::string {FLAN_SOURCE_DIR}
		+ "/data/sputtering/D_on_W.h5"};

	/// A temporary directory for a test background, removed at the end
	std::filesystem::path temp_dir(const std::string& name)
	{
		return std::filesystem::temp_directory_path() / ("flan_test_"
			+ name + "_" + std::to_string(getpid()));
	}

	Options::Options tungsten_options()
	{
		Options::Options opts {};
		opts.set_imp_atom_num(74);
		opts.set_imp_num(1000);
		return opts;
	}

	/// Relative difference, scaled by the larger magnitude (or 1e-300)
	double rel(const double a, const double b)
	{
		return std::abs(a - b) / std::max({std::abs(a), std::abs(b),
			1e-300});
	}

	/**
	* Fills n dead slots on the GPU from the source and compares each with
	* the CPU's birth from slot i's generator.
	*/
	void births_match_cpu(Background::Background& bkg,
		const Sputtering::Source& source, const int n)
	{
		auto opts {tungsten_options()};
		opts.set_imp_num(n);

		Slots::Slots slots(n);
		slots.set_mass(183.84 * 1.66053906660e-27);
		slots.set_Z(74);
		for (int i {}; i < n; ++i)
		{
			slots.set_state(i, 1);
			slots.set_prompt_by(i, 5.0);	// overwritten by the fill
			slots.set_generation(i, 3);
		}

		Background::BackgroundDevice bkg_d {bkg.to_device()};
		Slots::SlotsDevice slots_d {slots.to_device()};
		pcg32* rngs_d {init_rngs_cuda(slots_d, seed)};
		Options::OptionsDevice* opts_d {opts.to_device()};
		Sputtering::SourceDevice source_d {source.to_device(0, n)};
		int rem {n};
		int alive {};
		Slots::fill_slots_gpu(slots_d, rem, alive, rngs_d, opts_d, bkg_d,
			&source_d);
		ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);
		EXPECT_EQ(rem, 0);
		EXPECT_EQ(alive, n);
		const Slots::Slots gpu {slots.to_host(slots_d)};
		Sputtering::free_source(source_d);
		EXPECT_EQ(source_d.device_id, -1);
		Options::free_opts(opts_d);
		free_rngs(rngs_d, slots_d.device_id);
		Slots::free_slots(slots_d);
		Background::free_bkg(bkg_d);

		int bad {};
		for (int i {}; i < n; ++i)
		{
			pcg32 rng(seed + i, i);
			const Sputtering::Birth b {source.sample_birth(rng)};
			const double speed {std::hypot(b.vX, b.vY, b.vZ)};
			const bool ok {gpu.state()[i] == 0 && gpu.zone()[i] == b.zone
				&& gpu.q()[i] == 0 && gpu.prompt_by()[i] == -1.0
				&& gpu.generation()[i] == 0
				&& rel(gpu.weight()[i], source.source_rate() / n) < 1e-15
				&& rel(gpu.t()[i], b.t) < 1e-12
				&& rel(gpu.x()[i], b.x) < 1e-12
				&& rel(gpu.y()[i], b.y) < 1e-12
				&& std::abs(gpu.z()[i] - b.z) < 1e-12
				&& std::abs(gpu.vX()[i] - b.vX) < 1e-9 * speed
				&& std::abs(gpu.vY()[i] - b.vY) < 1e-9 * speed
				&& std::abs(gpu.vZ()[i] - b.vZ) < 1e-9 * speed};
			if (!ok && bad++ < 5)
			{
				ADD_FAILURE() << "birth " << i << ": GPU zone "
					<< gpu.zone()[i] << " (" << gpu.x()[i] << ", "
					<< gpu.y()[i] << ", " << gpu.z()[i] << ") v ("
					<< gpu.vX()[i] << ", " << gpu.vY()[i] << ", "
					<< gpu.vZ()[i] << "), CPU zone " << b.zone << " ("
					<< b.x << ", " << b.y << ", " << b.z << ") v (" << b.vX
					<< ", " << b.vY << ", " << b.vZ << ")";
			}
		}
		EXPECT_EQ(bad, 0);
	}
}

// HERMES targets, impact energy from the formula (a Gamma(2, Ti) thermal
// part, picked by its yield)
TEST(SputteringCUDA, BirthsMatchCPUOnHermesTargets)
{
	if (!std::filesystem::exists(d_on_w)) GTEST_SKIP() << d_on_w;
	const auto dir {temp_dir("sput_births_hermes")};
	Background::Background bkg {MultiZone::hermes_background(dir)};
	std::filesystem::remove_all(dir);

	const Sputtering::Surface surface(bkg);
	ASSERT_TRUE(surface.from_wall_faces());
	const Sputtering::YieldTable table(d_on_w);
	const Sputtering::Source source(surface, bkg, table, tungsten_options());
	ASSERT_GT(source.source_rate(), 0.0);
	births_match_cpu(bkg, source, 4000);
}

// SOLEDGE3X staircase wall, impact energies and angles from test ions
// pushed through the sheath (the impact samples)
TEST(SputteringCUDA, BirthsMatchCPUOnS3XStaircaseWithSheath)
{
	if (!std::filesystem::exists(d_on_w)) GTEST_SKIP() << d_on_w;
	const auto dir {temp_dir("sput_births_s3x")};
	Background::Background bkg {MultiZone::s3x_background(dir)};
	std::filesystem::remove_all(dir);

	const Sputtering::Surface surface(bkg);
	ASSERT_FALSE(surface.from_wall_faces());
	const Sputtering::YieldTable table(d_on_w);
	const Sputtering::Source source(surface, bkg, table, tungsten_options(),
		Sheath::Model::brooks);
	ASSERT_GT(source.source_rate(), 0.0);
	births_match_cpu(bkg, source, 4000);
}

// One ionization step with the charge limited to 2: the charges and the
// prompt windows started on the first ionization match the CPU's, nothing
// ionizes past 2, and some neutrals start a window
TEST(SputteringCUDA, IonizationPromptWindowAndChargeLimit)
{
	using namespace OpenADASTest;
	const auto root {write_fake_adas_dir("sput_cuda",
		[](int, double, double) {return 1e-14;},
		[](int, double, double) {return 3e-15;})};
	OpenADAS::OpenADAS oa_ioniz(root.string(), fake_year, fake_Z, "scd");
	OpenADAS::OpenADAS oa_recomb(root.string(), fake_year, fake_Z, "acd");

	Background::Background bkg {make_uniform_bkg(1e19, 10.0)};
	constexpr int per_charge {500};
	Slots::Slots slots {make_slots(bkg, per_charge)};
	slots.set_mass(6.94 * 1.66053906660e-27);
	pcg32 setup(3, 3);
	for (int i {}; i < slots.N(); ++i)
	{
		slots.set_pt(i, 1e-6 * setup.next_double());
		// A neutral that has been ionized before keeps its window
		if (i % 7 == 0) slots.set_prompt_by(i, 2e-6);
	}
	constexpr double dt {5e-6};		// ionization probability 0.5
	constexpr int max_charge {2};

	// CPU, one at a time with slot i's generator
	Slots::Slots cpu {slots};
	int iw {};
	int rw {};
	for (int i {}; i < slots.N(); ++i)
	{
		Slots::Slots one(1);
		one.set_mass(slots.mass());
		one.set_Z(slots.Z());
		one.set_state(0, 0);
		one.set_q(0, slots.q()[i]);
		one.set_t(0, slots.t()[i]);
		one.set_x(0, slots.x()[i]);
		one.set_y(0, slots.y()[i]);
		one.set_z(0, slots.z()[i]);
		one.set_pt(0, slots.pt()[i]);
		one.set_prompt_by(0, slots.prompt_by()[i]);
		std::vector<pcg32> rngs(omp_get_max_threads(), pcg32(seed + i, i));
		OpenADAS::ioniz_recomb(one, bkg, oa_ioniz, oa_recomb, dt, iw, rw,
			rngs, max_charge);
		cpu.set_q(i, one.q()[0]);
		cpu.set_prompt_by(i, one.prompt_by()[0]);
	}

	// GPU
	OpenADAS::OpenADASDevice oa_ioniz_d {oa_ioniz.to_device(0)};
	OpenADAS::OpenADASDevice oa_recomb_d {oa_recomb.to_device(0)};
	Background::BackgroundDevice bkg_d {bkg.to_device()};
	Slots::SlotsDevice slots_d {slots.to_device()};
	pcg32* rngs_d {init_rngs_cuda(slots_d, seed)};
	OpenADAS::IonizRecombWarningsDevice warnings_d {
		OpenADAS::alloc_warnings(0)};
	OpenADAS::ioniz_recomb_gpu(slots_d, bkg_d, oa_ioniz_d, oa_recomb_d, dt,
		warnings_d, rngs_d, max_charge);
	ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);
	const Slots::Slots gpu {slots.to_host(slots_d)};
	OpenADAS::free_warnings(warnings_d);
	free_rngs(rngs_d, slots_d.device_id);
	Slots::free_slots(slots_d);
	Background::free_bkg(bkg_d);
	OpenADAS::free_oa(oa_ioniz_d, 0);
	OpenADAS::free_oa(oa_recomb_d, 0);
	std::filesystem::remove_all(root);

	int started {};
	int bad {};
	for (int i {}; i < slots.N(); ++i)
	{
		EXPECT_LE(cpu.q()[i], std::max(slots.q()[i], max_charge));
		const bool ok {gpu.q()[i] == cpu.q()[i]
			&& rel(gpu.prompt_by()[i], cpu.prompt_by()[i]) < 1e-12};
		if (!ok && bad++ < 5)
		{
			ADD_FAILURE() << "particle " << i << ": GPU q " << gpu.q()[i]
				<< " prompt_by " << gpu.prompt_by()[i] << ", CPU q "
				<< cpu.q()[i] << " prompt_by " << cpu.prompt_by()[i];
		}
		if (slots.prompt_by()[i] < 0.0 && cpu.prompt_by()[i] >= 0.0)
		{
			++started;
			EXPECT_GT(cpu.prompt_by()[i], slots.pt()[i]);
		}
	}
	EXPECT_EQ(bad, 0);
	EXPECT_GT(started, per_charge / 4);

	// Nothing at the limit ionized
	for (int i {}; i < slots.N(); ++i)
	{
		if (slots.q()[i] == max_charge) EXPECT_LE(gpu.q()[i], max_charge);
	}
}

namespace
{
	/// Particles past the faces of their zones (as in
	/// test_boundary_zones_cuda), with random times since birth and
	/// prompt windows, so every tally gets something
	Slots::Slots landing_particles(const Background::Background& bkg,
		const int per_zone)
	{
		Slots::Slots slots {MultiZone::spread_particles(bkg, per_zone, 17)};
		pcg32 rng(8, 2);
		const auto& zones {bkg.get_zones()};
		for (int i {}; i < slots.N(); ++i)
		{
			const int zone {slots.zone()[i]};
			const double pos[3] {slots.x()[i], slots.y()[i], slots.z()[i]};
			double moved[3] {};
			for (int d {}; d < 3; ++d)
			{
				// A step of up to 30% of the zone either way from where it
				// was, so the slot's cell is the one it came from
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
			slots.set_weight(i, 1.0 + rng.next_double());
			if (i % 5 == 0) slots.set_q(i, 0);
		}
		return slots;
	}

	/**
	* Runs the boundary check with deposition on both and compares the
	* tallies.
	*
	* @return The CPU's tallies, for checks of what was tallied.
	*/
	Sputtering::Deposition deposition_matches_cpu(
		Background::Background& bkg, const Sputtering::Surface& surface,
		const int per_zone)
	{
		const Options::Options opts {};
		const Slots::Slots slots {landing_particles(bkg, per_zone)};

		// CPU, one at a time with slot i's generator
		Sputtering::Deposition cpu(surface, 74, bkg.get_phi_sign());
		Boundary::Fates fates {};
		for (int i {}; i < slots.N(); ++i)
		{
			Slots::Slots one(1);
			one.set_mass(slots.mass());
			one.set_Z(slots.Z());
			one.set_state(0, 0);
			one.set_zone(0, slots.zone()[i]);
			one.set_t(0, slots.t()[i]);
			one.set_x(0, slots.x()[i]);
			one.set_y(0, slots.y()[i]);
			one.set_z(0, slots.z()[i]);
			one.set_xidx(0, slots.xidx()[i]);
			one.set_yidx(0, slots.yidx()[i]);
			one.set_zidx(0, slots.zidx()[i]);
			one.set_vX(0, slots.vX()[i]);
			one.set_vY(0, slots.vY()[i]);
			one.set_vZ(0, slots.vZ()[i]);
			one.set_q(0, slots.q()[i]);
			one.set_pt(0, slots.pt()[i]);
			one.set_prompt_by(0, slots.prompt_by()[i]);
			one.set_weight(0, slots.weight()[i]);
			std::vector<pcg32> rngs(omp_get_max_threads(),
				pcg32(seed + i, i));
			Boundary::check_bounds_cpu(one, bkg, opts, rngs, fates, &cpu);
		}

		// GPU, tallied on the device and added to empty host tallies
		Sputtering::Deposition gpu(surface, 74, bkg.get_phi_sign());
		Background::BackgroundDevice bkg_d {bkg.to_device()};
		Slots::Slots copy {slots};
		Slots::SlotsDevice slots_d {copy.to_device()};
		pcg32* rngs_d {init_rngs_cuda(slots_d, seed)};
		Sputtering::DepositionDevice dep_d {gpu.to_device(0)};
		Boundary::check_bounds_gpu(slots_d, bkg_d,
			Boundary::make_params(bkg, opts, bkg_d.chi), rngs_d, nullptr,
			&dep_d);
		EXPECT_EQ(cudaDeviceSynchronize(), cudaSuccess);
		gpu.add_from_device(dep_d);
		Sputtering::free_deposition(dep_d);
		free_rngs(rngs_d, slots_d.device_id);
		Slots::free_slots(slots_d);
		Background::free_bkg(bkg_d);

		auto same = [](const std::vector<double>& g,
			const std::vector<double>& c, const char* name)
		{
			ASSERT_EQ(g.size(), c.size()) << name;
			int bad {};
			for (std::size_t n {}; n < c.size(); ++n)
			{
				if (!(rel(g[n], c[n]) < 1e-12) && bad++ < 3)
				{
					ADD_FAILURE() << name << "[" << n << "]: GPU " << g[n]
						<< ", CPU " << c[n];
				}
			}
		};
		same(gpu.redeposited(), cpu.redeposited(), "redeposited");
		same(gpu.prompt(), cpu.prompt(), "prompt");
		same(gpu.element_energy(), cpu.element_energy(), "element_energy");
		same(gpu.lifetime(), cpu.lifetime(), "lifetime");
		same(gpu.charge(), cpu.charge(), "charge");
		same(gpu.prompt_charge(), cpu.prompt_charge(), "prompt_charge");
		same(gpu.impact_energy(), cpu.impact_energy(), "impact_energy");
		same(gpu.impact_angle(), cpu.impact_angle(), "impact_angle");
		EXPECT_LT(rel(gpu.unassigned(), cpu.unassigned()), 1e-12);
		EXPECT_LT(rel(gpu.wall_face(), cpu.wall_face()), 1e-12);
		EXPECT_LT(rel(gpu.mean_lifetime(), cpu.mean_lifetime()), 1e-12);
		EXPECT_LT(rel(gpu.mean_energy(), cpu.mean_energy()), 1e-12);
		EXPECT_LT(rel(gpu.mean_angle(), cpu.mean_angle()), 1e-12);
		return cpu;
	}

	double sum(const std::vector<double>& v)
	{
		double s {};
		for (const double x : v) s += x;
		return s;
	}
}

// HERMES targets: atoms land on the element at the cell of the face they
// cross
TEST(SputteringCUDA, DepositionMatchesCPUAtHermesTargets)
{
	const auto dir {temp_dir("sput_dep_hermes")};
	Background::Background bkg {MultiZone::hermes_background(dir)};
	std::filesystem::remove_all(dir);
	const Sputtering::Surface surface(bkg);
	ASSERT_TRUE(surface.from_wall_faces());

	const auto cpu {deposition_matches_cpu(bkg, surface, 2000)};
	EXPECT_GT(sum(cpu.redeposited()), 0.0);
	EXPECT_GT(sum(cpu.prompt()), 0.0);
	EXPECT_GT(sum(cpu.impact_angle()), 0.0);
}

// SOLEDGE3X: atoms absorbed in wall cells go on the staircase element
// between the wall cell and the cell they came from; at the wall faces
// (inside the material) they're only counted there
TEST(SputteringCUDA, DepositionMatchesCPUInS3XWallCells)
{
	const auto dir {temp_dir("sput_dep_s3x")};
	Background::Background bkg {MultiZone::s3x_background(dir)};
	std::filesystem::remove_all(dir);
	const Sputtering::Surface surface(bkg);
	ASSERT_FALSE(surface.from_wall_faces());

	const auto cpu {deposition_matches_cpu(bkg, surface, 2000)};
	EXPECT_GT(sum(cpu.redeposited()), 0.0);
	EXPECT_GT(sum(cpu.prompt()), 0.0);
	EXPECT_GT(cpu.wall_face() + cpu.unassigned(), 0.0);
}
