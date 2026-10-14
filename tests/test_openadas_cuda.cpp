#include <array>
#include <gtest/gtest.h>
#include <vector>

#include <cuda_runtime.h>

#include "background.h"
#include "background_device.h"
#include "openadas.h"
#include "openadas_device.h"
#include "pcg32.h"
#include "slots.h"
#include "slots_device.h"

#include "init_rngs.cuh"
#include "openadas.cuh"

#include "openadas_test_utils.h"

using namespace OpenADASTest;

namespace
{
	struct StepResult
	{
		int ioniz_warnings {};
		int recomb_warnings {};
	};

	// Loads the fake ADAS files in root, copies everything to the GPU, takes
	// nsteps ionization/recombination steps with ioniz_recomb_gpu and copies
	// the slots back to the host. The warning counters are allocated once
	// and accumulate over all the steps, the same as in a simulation, and
	// are read back at the end.
	StepResult step_gpu(const fs::path& root,
		Background::Background& bkg, Slots::Slots& slots,
		const double dt, uint64_t seed = 4, int nsteps = 1)
	{
		OpenADAS::OpenADAS oa_ioniz(root.string(), fake_year, fake_Z, "scd");
		OpenADAS::OpenADAS oa_recomb(root.string(), fake_year, fake_Z, "acd");
		OpenADAS::OpenADASDevice oa_ioniz_d {oa_ioniz.to_device(0)};
		OpenADAS::OpenADASDevice oa_recomb_d {oa_recomb.to_device(0)};

		Background::BackgroundDevice bkg_d {bkg.to_device()};
		Slots::SlotsDevice slots_d {slots.to_device()};
		pcg32* rngs_d {init_rngs_cuda(slots_d, seed)};
		OpenADAS::IonizRecombWarningsDevice warnings_d {
			OpenADAS::alloc_warnings(0)};

		for (int n {}; n < nsteps; ++n)
		{
			OpenADAS::ioniz_recomb_gpu(slots_d, bkg_d, oa_ioniz_d,
				oa_recomb_d, dt, warnings_d, rngs_d);
		}

		StepResult res {};
		OpenADAS::add_warnings_to_host(warnings_d, res.ioniz_warnings,
			res.recomb_warnings);

		// The kernel launches are asynchronous now, so check here that they
		// all ran without error.
		EXPECT_EQ(cudaDeviceSynchronize(), cudaSuccess);

		// Copy back to host and free memory on the GPU
		slots = slots.to_host(slots_d);
		Slots::free_slots(slots_d);
		Background::free_bkg(bkg_d);
		cudaFree(rngs_d);
		OpenADAS::free_warnings(warnings_d);
		OpenADAS::free_oa(oa_ioniz_d, 0);
		OpenADAS::free_oa(oa_recomb_d, 0);

		return res;
	}
}


// -------------------------------------------------------------------------
// Charge-state update (deterministic cases)
// -------------------------------------------------------------------------

// With certain ionization and no recombination, every particle that can
// ionize goes up by one. The fully stripped particles can't ionize.
TEST(IonizRecomb, IonizationOnly_IncreasesCharge)
{
	const fs::path root {write_fake_adas_dir("ioniz_only_cuda", huge_rate,
		tiny_rate)};
	Background::Background bkg {make_uniform_bkg(1e19, 10.0)};

	constexpr int per_charge {16};
	Slots::Slots slots {make_slots(bkg, per_charge)};
	StepResult res {step_gpu(root, bkg, slots, 1e-6)};

	for (int i {}; i < slots.N(); ++i)
	{
		const int q0 {i / per_charge};
		EXPECT_EQ(slots.q()[i], std::min(q0 + 1, fake_Z)) << "q0 = " << q0;
	}

	// Probability > 1 for every particle that can ionize
	EXPECT_EQ(res.ioniz_warnings, fake_Z * per_charge);
	EXPECT_EQ(res.recomb_warnings, 0);
}


// With certain recombination and no ionization, every particle that can
// recombine goes down by one. Neutrals can't recombine.
TEST(IonizRecomb, RecombinationOnly_DecreasesCharge)
{
	const fs::path root {write_fake_adas_dir("recomb_only_cuda", tiny_rate,
		huge_rate)};
	Background::Background bkg {make_uniform_bkg(1e19, 10.0)};

	constexpr int per_charge {16};
	Slots::Slots slots {make_slots(bkg, per_charge)};
	StepResult res {step_gpu(root, bkg, slots, 1e-6)};

	for (int i {}; i < slots.N(); ++i)
	{
		const int q0 {i / per_charge};
		EXPECT_EQ(slots.q()[i], std::max(q0 - 1, 0)) << "q0 = " << q0;
	}

	EXPECT_EQ(res.ioniz_warnings, 0);
	EXPECT_EQ(res.recomb_warnings, fake_Z * per_charge);
}


// With both certain, they cancel for partially ionized particles. Neutrals
// can only ionize and fully stripped particles can only recombine.
TEST(IonizRecomb, Cancellation_NoNetChange)
{
	const fs::path root {write_fake_adas_dir("cancel_cuda", huge_rate,
		huge_rate)};
	Background::Background bkg {make_uniform_bkg(1e19, 10.0)};

	constexpr int per_charge {16};
	Slots::Slots slots {make_slots(bkg, per_charge)};
	StepResult res {step_gpu(root, bkg, slots, 1e-6)};

	for (int i {}; i < slots.N(); ++i)
	{
		const int q0 {i / per_charge};
		int expected {q0};
		if (q0 == 0) expected = 1;
		else if (q0 == fake_Z) expected = fake_Z - 1;
		EXPECT_EQ(slots.q()[i], expected) << "q0 = " << q0;
	}

	EXPECT_EQ(res.ioniz_warnings, fake_Z * per_charge);
	EXPECT_EQ(res.recomb_warnings, fake_Z * per_charge);
}


// The warning counters live on the GPU for the whole run, so they must keep
// adding up across steps rather than reset each step. With certain
// ionization, each step every particle below fully stripped ionizes and
// warns: 3 charge groups warn on step 1, 2 on step 2 and 1 on step 3, after
// which everything is fully stripped and stays that way.
TEST(IonizRecomb, Warnings_AccumulateAcrossSteps)
{
	const fs::path root {write_fake_adas_dir("accumulate_cuda", huge_rate,
		tiny_rate)};
	Background::Background bkg {make_uniform_bkg(1e19, 10.0)};

	constexpr int per_charge {16};
	Slots::Slots slots {make_slots(bkg, per_charge)};
	StepResult res {step_gpu(root, bkg, slots, 1e-6, 4, 5)};

	for (int i {}; i < slots.N(); ++i)
		EXPECT_EQ(slots.q()[i], fake_Z);

	EXPECT_EQ(res.ioniz_warnings, (3 + 2 + 1) * per_charge);
	EXPECT_EQ(res.recomb_warnings, 0);
}


// Dead particles are skipped entirely, even with certain ionization and
// recombination.
TEST(IonizRecomb, DeadParticles_Unchanged)
{
	const fs::path root {write_fake_adas_dir("dead_cuda", huge_rate,
		huge_rate)};
	Background::Background bkg {make_uniform_bkg(1e19, 10.0)};

	constexpr int per_charge {16};
	Slots::Slots slots {make_slots(bkg, per_charge)};
	for (int i {}; i < slots.N(); ++i) slots.set_state(i, 1);
	StepResult res {step_gpu(root, bkg, slots, 1e-6)};

	for (int i {}; i < slots.N(); ++i)
		EXPECT_EQ(slots.q()[i], i / per_charge);

	EXPECT_EQ(res.ioniz_warnings, 0);
	EXPECT_EQ(res.recomb_warnings, 0);
}


// -------------------------------------------------------------------------
// Charge-state update (statistical)
// -------------------------------------------------------------------------

// Many particles in each charge state take one step at a given (ne, Te).
// The fraction that ionize and recombine must match the analytic rates.
// This is the GPU check that the rate table is read as (charge, Te, ne),
// with q for ionization and q-1 for recombination, from tables that use
// charge_low = 1 and are non-square, like the real ADF11 files. The two
// plasma conditions are off-grid points where ionization dominates and
// where recombination is significant.
TEST(IonizRecomb, FractionsMatchRates)
{
	const fs::path root {write_fake_adas_dir("fractions_cuda", scd_rate,
		acd_rate)};

	// ne (m^-3), Te (eV), dt (s)
	const std::vector<std::array<double, 3>> cases {
		{4e19, 30.0, 1e-8}, {5e20, 3.0, 5e-10}};

	constexpr int per_charge {100000};
	uint64_t seed {11};
	for (auto [ne, te, dt] : cases)
	{
		Background::Background bkg {make_uniform_bkg(ne, te)};
		Slots::Slots slots {make_slots(bkg, per_charge)};
		StepResult res {step_gpu(root, bkg, slots, dt, seed++)};

		expect_fractions_match(slots, per_charge, ne, te, dt);
		EXPECT_EQ(res.ioniz_warnings, 0);
		EXPECT_EQ(res.recomb_warnings, 0);
	}
}
