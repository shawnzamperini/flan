#include <gtest/gtest.h>
#include "slots.h"


TEST(FillSlotsCPU, RevivesDeadParticles)
{
	// 10 slots, 5 dead particles and only 3 particles left to swap in
    Slots::Slots slots(10);
    for (int i = 0; i < 10; i++) slots.set_state(i, i % 2); // 5 dead
    int rem = 3;
	int alive {};
    Slots::fill_slots_cpu(slots, rem, alive);

	// 10 alive, then minus 5 dead, then 3 alive ones swapped in so we expect
	// 8 alive particles and 0 remaining ones
    EXPECT_EQ(alive, 8);
    EXPECT_EQ(rem, 0);
}
