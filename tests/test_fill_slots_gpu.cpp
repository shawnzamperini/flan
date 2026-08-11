#include <gtest/gtest.h>
#include "slots.h"
#include "slots_device.h"
#include "slots.cuh"

TEST(FillSlotsGPU, RevivesDeadParticles)
{
#ifdef USE_CUDA
std::cout << "CUDA enabled\n";
#else
std::cout << "CUDA disabled\n";
#endif

	std::cout << "Test starting...\n";
	// 10 slots, 5 dead particles and only 3 particles left to swap in
    Slots::Slots slots(10);
    for (int i = 0; i < 10; i++) slots.set_state(i, i % 2); // 5 dead
    int rem = 3;

	std::cout << "before\n";
	for (int i = 0; i < 10; i++) std::cout << i << " : " << slots.state()[i] << '\n';

	// Copy to device
	Slots::SlotsDevice slots_d {slots.to_device()};

	// Call wrapper that calls kernel to fill slots on the GPU
	Slots::fill_slots_gpu(slots_d, rem);

	// Copy back to host, copies slots_d into slots. Free memory on GPU.
	slots = slots.to_host(slots_d);
	Slots::free_slots(slots_d);

	// Count number of alive particles in slots
    int alive = 0;
    for (int i = 0; i < 10; i++)
        if (slots.state()[i] == 0) alive++;

	std::cout << "after\n";
	for (int i = 0; i < 10; i++) std::cout << i << " : " << slots.state()[i] << '\n';

	// 10 alive, then minus 5 dead, then 3 alive ones swapped in so we expect
	// 8 alive particles and 0 remaining ones
    EXPECT_EQ(alive, 8);
    EXPECT_EQ(rem, 0);
}
