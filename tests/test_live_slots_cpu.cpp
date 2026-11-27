/*
* How the GPU kernels' threads are given the live slots
* (Slots::per_warp_for, Slots::list_position, Slots::num_blocks): every
* entry of the live list goes to exactly one thread of the launch, warps
* get at most per_warp of them, and the launch fills the GPU's resident
* warps while there are fewer live slots than that. Only index arithmetic,
* so it runs without a GPU.
*/

#include <gtest/gtest.h>
#include <vector>

#include "slots_device.h"

namespace
{
	/// A SlotsDevice with just the fields num_blocks reads
	Slots::SlotsDevice with_list(const int num_live, const int per_warp)
	{
		Slots::SlotsDevice s {};
		s.num_live = num_live;
		s.per_warp = per_warp;
		return s;
	}
}

// Few live slots get a warp each, many are packed 32 to a warp
TEST(LiveSlotsCPU, PerWarp)
{
	EXPECT_EQ(Slots::per_warp_for(100, 0), Slots::warp_size);
	EXPECT_EQ(Slots::per_warp_for(0, 864), 1);
	EXPECT_EQ(Slots::per_warp_for(1, 864), 1);
	EXPECT_EQ(Slots::per_warp_for(864, 864), 1);
	EXPECT_EQ(Slots::per_warp_for(865, 864), 2);
	EXPECT_EQ(Slots::per_warp_for(4000, 864), 5);
	EXPECT_EQ(Slots::per_warp_for(32 * 864, 864), 32);
	EXPECT_EQ(Slots::per_warp_for(1000000, 864), 32);
}

// Every list entry goes to exactly one thread, with no more than
// per_warp to a warp, for both block sizes the launchers use
TEST(LiveSlotsCPU, EachEntryOnce)
{
	for (const int block_size : {128, 256})
	for (const int num_live : {1, 31, 32, 33, 500, 864, 4000, 100000})
	{
		const int per_warp {Slots::per_warp_for(num_live, 864)};
		const int threads {Slots::num_blocks(with_list(num_live, per_warp),
			block_size) * block_size};
		std::vector<int> times (num_live);
		std::vector<int> in_warp ((threads + 31) / 32);
		for (int k {}; k < threads; ++k)
		{
			const int j {Slots::list_position(k, per_warp, num_live)};
			if (j < 0) continue;
			ASSERT_LT(j, num_live);
			++times[j];
			++in_warp[k / Slots::warp_size];
		}
		for (int j {}; j < num_live; ++j)
			ASSERT_EQ(times[j], 1) << "entry " << j << " of " << num_live;
		for (const int n : in_warp) EXPECT_LE(n, per_warp);

		// No more warps with work than the GPU holds until they're full,
		// and no more than a block's worth of threads left over
		const int busy_warps {(num_live + per_warp - 1) / per_warp};
		if (per_warp < Slots::warp_size) EXPECT_LE(busy_warps, 864);
		EXPECT_LT(threads - block_size,
			(num_live + per_warp - 1) / per_warp * Slots::warp_size);
	}
}

// Packed fully (per_warp 32), thread k takes entry k, as before spreading
TEST(LiveSlotsCPU, PackedIsInOrder)
{
	for (int k {}; k < 1000; ++k)
		EXPECT_EQ(Slots::list_position(k, Slots::warp_size, 900),
			k < 900 ? k : -1);
}
