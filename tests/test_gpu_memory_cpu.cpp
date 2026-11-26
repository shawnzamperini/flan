/*
* The check that stops a GPU run whose background doesn't fit in the GPU's
* memory (ImpurityTransport::gpu_memory_shortfall): a run that fits with
* room to spare gives no message, and one that doesn't says how much it
* needs and how much the GPU has. Runs without a GPU, since it only does
* arithmetic on the sizes.
*/

#include <cstddef>
#include <gtest/gtest.h>
#include <string>

#include "impurity_transport.h"

namespace
{
	constexpr std::size_t GiB {std::size_t {1} << 30};
	constexpr std::size_t MiB {std::size_t {1} << 20};
}

// Fits: 20 GiB on a GPU with 39 GiB free
TEST(GpuMemory, FitsWithRoom)
{
	ImpurityTransport::GpuMemoryNeed need {};
	need.background = 12 * GiB;
	need.statistics = 7 * GiB;
	need.slots = 1 * GiB;
	EXPECT_EQ(ImpurityTransport::gpu_memory_shortfall(need, 39 * GiB,
		40 * GiB, 0), "");
}

// Doesn't fit: the message has every size and the GPU
TEST(GpuMemory, TooBig)
{
	ImpurityTransport::GpuMemoryNeed need {};
	need.background = 45 * GiB;
	need.statistics = 25 * GiB;
	need.slots = GiB / 4;
	const std::string msg {ImpurityTransport::gpu_memory_shortfall(need,
		39 * GiB, 40 * GiB, 2)};
	ASSERT_FALSE(msg.empty());
	EXPECT_NE(msg.find("GPU 2"), std::string::npos) << msg;
	EXPECT_NE(msg.find("70.25 GiB"), std::string::npos) << msg;
	EXPECT_NE(msg.find("background 45.00 GiB"), std::string::npos) << msg;
	EXPECT_NE(msg.find("statistics 25.00 GiB"), std::string::npos) << msg;
	EXPECT_NE(msg.find("39.00 GiB free of 40.00 GiB"), std::string::npos)
		<< msg;
}

// Fits on its own but not with the room kept for the smaller copies: a
// tenth of the free memory, at least 256 MiB
TEST(GpuMemory, KeepsRoomToSpare)
{
	ImpurityTransport::GpuMemoryNeed need {};
	need.background = 36 * GiB;
	EXPECT_FALSE(ImpurityTransport::gpu_memory_shortfall(need, 39 * GiB,
		40 * GiB, 0).empty());
	need.background = 35 * GiB;
	EXPECT_TRUE(ImpurityTransport::gpu_memory_shortfall(need, 39 * GiB,
		40 * GiB, 0).empty());

	// A small GPU: at least 256 MiB to spare
	need.background = 900 * MiB;
	EXPECT_FALSE(ImpurityTransport::gpu_memory_shortfall(need, GiB, GiB,
		0).empty());
	need.background = 700 * MiB;
	EXPECT_TRUE(ImpurityTransport::gpu_memory_shortfall(need, GiB, GiB,
		0).empty());
}
