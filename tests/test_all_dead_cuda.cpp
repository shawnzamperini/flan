#include <gtest/gtest.h>
#include "slots.h"
#include "slots_device.h"
#include "slots.cuh"


TEST(AllDeadCUDA, AllDead)
{
	// 10 slots, all dead
    Slots::Slots slots(10);
    for (int i = 0; i < 10; i++) slots.set_state(i, 1); // all dead

	// Copy to device
	Slots::SlotsDevice slots_d {slots.to_device()};

	// Call wrapper that checks if all particles are dead
	bool all_dead {Slots::all_dead_gpu(slots_d)};

	// Free memory on GPU
	Slots::free_slots(slots_d);

	// Should be true
	EXPECT_TRUE(all_dead);
}

TEST(AllDeadCUDA, NotAllDead)
{
	// 10 slots, all dead
    Slots::Slots slots(10);
    for (int i = 0; i < 10; i++) slots.set_state(i, 1); // all dead

	// Set first slot to alive
	slots.set_state(0, 0);

	// Copy to device
	Slots::SlotsDevice slots_d {slots.to_device()};

	// Call wrapper that checks if all particles are dead
	bool all_dead {Slots::all_dead_gpu(slots_d)};

	// Free memory on GPU
	Slots::free_slots(slots_d);

	// Should be false
	EXPECT_FALSE(all_dead);
}
