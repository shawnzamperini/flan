#include <gtest/gtest.h>

#include "options.h"
#include "options_device.h"
#include "pcg32.h"
#include "slots.h"
#include "slots_device.h"

#include "init_rngs.cuh"
#include "slots.cuh"

TEST(FillSlotsCUDA, RevivesDeadParticles)
{
	// 10 slots, 5 dead particles and only 3 particles left to swap in
    Slots::Slots slots(10);
    for (int i = 0; i < 10; i++) slots.set_state(i, i % 2); // 5 dead
    int rem = 3;

	// Copy to device
	Slots::SlotsDevice slots_d {slots.to_device()};

	// Create PCG32 RNGs with just base_seed = 1. Unimportant for this test.
	pcg32* rngs_d {init_rngs_cuda(slots_d, 1)};

	// Create dummy options to pass in. Options not important for this test.
	Options::Options opts {};
	Options::OptionsDevice* opts_d {opts.to_device()};

	// Call wrapper that calls kernel to fill slots on the GPU
	int alive {};
	Slots::fill_slots_gpu(slots_d, rem, alive, rngs_d, opts_d);

	// Copy back to host, copies slots_d into slots. Free memory on GPU.
	slots = slots.to_host(slots_d);
	Slots::free_slots(slots_d);
	Options::free_opts(opts_d);

	// 10 alive, then minus 5 dead, then 3 alive ones swapped in so we expect
	// 8 alive particles and 0 remaining ones
    EXPECT_EQ(alive, 8);
    EXPECT_EQ(rem, 0);
}
