/*
* The list of live slots that fill_slots_gpu makes for the kernels
* (SlotsDevice::live_indices, Slots::slot_of_thread): it holds exactly the
* live slots in increasing order after a fill, and running the kernels
* over it gives the same slots, bit for bit, as running them over every
* slot. Each slot keeps its own generator, so which thread handles it
* mustn't matter.
*
* The slots are set up as in a run that has stopped refilling: every slot
* filled, most of them then killed at scattered places, and a few refilled.
*/

#include <cstdint>
#include <cuda_runtime.h>
#include <gtest/gtest.h>
#include <vector>

#include "background.h"
#include "background_device.h"
#include "impurity_stats.h"
#include "impurity_stats_device.h"
#include "options.h"
#include "options_device.h"
#include "pcg32.h"
#include "read_bkg.h"
#include "slots.h"
#include "slots_device.h"

#include "boris.cuh"
#include "collisions.cuh"
#include "impurity_transport.cuh"
#include "init_rngs.cuh"
#include "slots.cuh"

namespace
{
	constexpr int num_slots {3000};
	constexpr int refill {200};

	/// Whether slot i is killed after the first fill: about 70%, in runs of
	/// different lengths so warps have a few live slots each
	bool killed(const int i)
	{
		const std::uint32_t h {static_cast<std::uint32_t>(i) * 2654435761u};
		return (h >> 16) % 10 < 7;
	}

	/// Slots, generators and background on the GPU, set up as described
	/// at the top
	struct LiveSetup
	{
		Options::Options opts {};
		Background::Background bkg {};
		Background::BackgroundDevice bkg_d {};
		Options::OptionsDevice* opts_d {nullptr};
		Slots::Slots slots {num_slots};
		Slots::SlotsDevice slots_d {};
		pcg32* rngs_d {nullptr};
		int alive {};
		int rem {};

		LiveSetup()
		{
			opts.set_bkg_source("slab");
			bkg = Background::read_bkg(opts);
			bkg_d = bkg.to_device();
			opts_d = opts.to_device();
			for (int i {}; i < num_slots; ++i) slots.set_state(i, 1);
			slots_d = slots.to_device();
			rngs_d = init_rngs_cuda(slots_d, 7);

			// Fill every slot
			rem = num_slots;
			Slots::fill_slots_gpu(slots_d, rem, alive, rngs_d, opts_d, bkg_d);

			// Kill most of them, scattered
			std::vector<int> state (num_slots);
			for (int i {}; i < num_slots; ++i) state[i] = killed(i) ? 1 : 0;
			cudaMemcpy(slots_d.state, state.data(), num_slots * sizeof(int),
				cudaMemcpyHostToDevice);

			// Refill a few
			rem = refill;
			Slots::fill_slots_gpu(slots_d, rem, alive, rngs_d, opts_d, bkg_d);
		}

		~LiveSetup()
		{
			free_rngs(rngs_d, slots_d.device_id);
			Slots::free_slots(slots_d);
			Options::free_opts(opts_d);
			Background::free_bkg(bkg_d, 0);
		}
	};

	/// Runs the push, both collision launches, the step and the cell search
	/// on the slots
	void run_kernels(Slots::SlotsDevice& slots_d, LiveSetup& s, pcg32* rngs_d)
	{
		Impurity::Statistics stats {s.bkg};
		ImpurityStats::StatisticsDevice stats_d {stats.to_device(0)};
		const double dt {1e-8};
		int* max_time_d {nullptr};
		cudaMalloc(&max_time_d, sizeof(int));
		cudaMemset(max_time_d, 0, sizeof(int));
		for (int step {}; step < 3; ++step)
		{
			for (const bool elec : {false, true})
			{
				Collisions::collision_gpu(slots_d, s.bkg_d, elec, dt, stats_d,
					rngs_d, s.bkg.get_elec_mass_amu(),
					s.bkg.get_ion_mass_amu(), s.opts_d, nullptr);
			}
			Boris::update_velocity_gpu(slots_d, s.bkg_d, dt);
			ImpurityTransport::step_gpu(slots_d, dt, s.opts_d, max_time_d);
			ImpurityTransport::find_containing_cell_gpu(slots_d, s.bkg_d);
		}
		cudaDeviceSynchronize();
		cudaFree(max_time_d);
		Impurity::free_stats(stats_d, 0);
	}

	template <typename T>
	void expect_same(const std::vector<T>& a, const std::vector<T>& b,
		const char* what)
	{
		ASSERT_EQ(a.size(), b.size()) << what;
		int differ {};
		for (std::size_t i {}; i < a.size(); ++i) differ += (a[i] != b[i]);
		EXPECT_EQ(differ, 0) << what << " differs in " << differ << " slots";
	}
}

// No list before the first fill: the kernels run over every slot
TEST(LiveSlotsCUDA, NoListBeforeFill)
{
	Slots::Slots slots(100);
	Slots::SlotsDevice slots_d {slots.to_device()};
	EXPECT_EQ(slots_d.live_indices, nullptr);
	EXPECT_EQ(slots_d.num_live, 100);
	Slots::free_slots(slots_d);
}

// After the fill, the list is exactly the live slots in increasing order,
// and its length is the host's count
TEST(LiveSlotsCUDA, ListsLiveSlotsInOrder)
{
	LiveSetup s {};
	ASSERT_NE(s.slots_d.live_indices, nullptr);
	EXPECT_EQ(s.slots_d.num_live, s.alive);
	EXPECT_EQ(s.rem, 0);

	int num_live {};
	cudaMemcpy(&num_live, s.slots_d.num_live_d, sizeof(int),
		cudaMemcpyDeviceToHost);
	EXPECT_EQ(num_live, s.alive);
	std::vector<int> live (num_live);
	cudaMemcpy(live.data(), s.slots_d.live_indices, num_live * sizeof(int),
		cudaMemcpyDeviceToHost);

	const Slots::Slots host {s.slots.to_host(s.slots_d)};
	std::vector<int> expected {};
	for (int i {}; i < num_slots; ++i)
		if (host.state()[i] == 0) expected.push_back(i);
	EXPECT_EQ(live, expected);

	// The setup gives a scattered list: most slots dead, some refilled
	EXPECT_GT(num_live, refill);
	EXPECT_LT(num_live, num_slots / 2);
}

// The kernels give the same slots over the list as over every slot
TEST(LiveSlotsCUDA, SameResultsWithAndWithoutList)
{
	LiveSetup s {};

	// A second copy of the same slots and generators, without a list
	const Slots::Slots start {s.slots.to_host(s.slots_d)};
	Slots::Slots copy {start};
	Slots::SlotsDevice all_d {copy.to_device()};
	ASSERT_EQ(all_d.live_indices, nullptr);
	pcg32* rngs_all {nullptr};
	cudaMalloc(&rngs_all, num_slots * sizeof(pcg32));
	cudaMemcpy(rngs_all, s.rngs_d, num_slots * sizeof(pcg32),
		cudaMemcpyDeviceToDevice);

	run_kernels(s.slots_d, s, s.rngs_d);
	run_kernels(all_d, s, rngs_all);

	const Slots::Slots a {s.slots.to_host(s.slots_d)};
	const Slots::Slots b {copy.to_host(all_d)};
	expect_same(a.x(), b.x(), "x");
	expect_same(a.y(), b.y(), "y");
	expect_same(a.z(), b.z(), "z");
	expect_same(a.vX(), b.vX(), "vX");
	expect_same(a.vY(), b.vY(), "vY");
	expect_same(a.vZ(), b.vZ(), "vZ");
	expect_same(a.t(), b.t(), "t");
	expect_same(a.state(), b.state(), "state");
	expect_same(a.xidx(), b.xidx(), "xidx");
	expect_same(a.zidx(), b.zidx(), "zidx");

	// And the kernels did something to the live slots
	int moved {};
	for (int i {}; i < num_slots; ++i)
		moved += (a.state()[i] == 0 && a.x()[i] != start.x()[i]);
	EXPECT_GT(moved, refill);

	cudaFree(rngs_all);
	Slots::free_slots(all_d);
}
