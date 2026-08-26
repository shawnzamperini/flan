#include <gtest/gtest.h>
#include <omp.h>
#include "pcg32.h"
#include "slots.h"


TEST(FillSlotsCPU, RevivesDeadParticles)
{
	// 10 slots, 5 dead particles and only 3 particles left to swap in
    Slots::Slots slots(10);
    for (int i = 0; i < 10; i++) slots.set_state(i, i % 2); // 5 dead
    int rem = 3;
	int alive {};
    Options::Options opts {};
    Background::Background bkg {};

	// Create RNGs
	std::vector<pcg32> rngs;
	int nthreads = omp_get_max_threads();
	rngs.reserve(nthreads);

	// Eventually make base_seed an input option
	constexpr uint64_t base_seed {4};
	constexpr uint64_t rank {0};
	for (int tid = 0; tid < nthreads; ++tid) 
	{
		uint64_t seed = base_seed + rank;  // unique per rank
		uint64_t stream = tid;             // unique per thread
		rngs.push_back(pcg32(seed, stream));
	}

    Slots::fill_slots_cpu(slots, rem, alive, bkg, opts, rngs);

	// 10 alive, then minus 5 dead, then 3 alive ones swapped in so we expect
	// 8 alive particles and 0 remaining ones
    EXPECT_EQ(alive, 8);
    EXPECT_EQ(rem, 0);
}
