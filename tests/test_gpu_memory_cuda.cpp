/*
* The GPU memory a background and its statistics take, as the run's check
* estimates it (Background::device_bytes, Statistics::device_bytes),
* against what copying them to the GPU takes (the drop in free memory from
* cudaMemGetInfo), on the multi-zone test backgrounds.
*
* cudaMemGetInfo counts memory in pages of 2 MiB, and these backgrounds
* are only a few MiB, so each is copied eight times and the drop compared
* with eight times the estimate, to 3% plus 2 MiB per copy. That catches
* a wrong count of frames, cells or ghost cells, not one array more or
* less; the estimate was checked to that on a large BOUT++ background
* (see the GPU plan, G8). Another process using the same GPU changes the
* free memory too, so the test runs on its own under ctest (RUN_SERIAL).
*/

#include <cstddef>
#include <cuda_runtime.h>
#include <filesystem>
#include <iostream>
#include <gtest/gtest.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "background.h"
#include "background_device.h"
#include "impurity_stats.h"
#include "impurity_stats_device.h"
#include "multizone_cases.h"

namespace
{
	constexpr int copies {8};

	std::size_t free_memory()
	{
		cudaDeviceSynchronize();
		std::size_t free_bytes {};
		std::size_t total_bytes {};
		cudaMemGetInfo(&free_bytes, &total_bytes);
		return free_bytes;
	}

	/// Copies the background and its statistics to GPU 0 several times
	/// and compares the memory taken with the estimate
	void estimate_matches(Background::Background& bkg)
	{
		cudaSetDevice(0);
		cudaFree(nullptr);	// start the context before measuring

		Impurity::Statistics stats {bkg};	// as follow_impurities makes it
		const std::size_t estimate {bkg.device_bytes()
			+ stats.device_bytes()};
		const std::size_t before {free_memory()};
		std::vector<Background::BackgroundDevice> bkgs;
		std::vector<ImpurityStats::StatisticsDevice> statss;
		for (int c {}; c < copies; ++c)
		{
			bkgs.push_back(bkg.to_device(0));
			statss.push_back(stats.to_device(0));
		}
		const std::size_t after {free_memory()};
		for (auto& b : bkgs) Background::free_bkg(b, 0);
		for (auto& s : statss) Impurity::free_stats(s, 0);

		const double used {static_cast<double>(before - after)};
		const double expected {static_cast<double>(copies * estimate)};
		std::cout << "estimate " << estimate / 1048576.0 << " MiB per "
			"copy, " << copies << " copies took " << used / 1048576.0
			<< " MiB\n";
		EXPECT_NEAR(used, expected, 0.03 * expected + copies * 2.0
			* 1048576.0);
	}

	std::filesystem::path temp_dir(const std::string& name)
	{
		return std::filesystem::temp_directory_path() / ("flan_test_"
			+ name + "_" + std::to_string(getpid()));
	}
}

TEST(GpuMemoryCUDA, EstimateMatchesHermes)
{
	const auto dir {temp_dir("gpu_mem_hermes")};
	Background::Background bkg {MultiZone::hermes_background(dir)};
	std::filesystem::remove_all(dir);
	estimate_matches(bkg);
}

TEST(GpuMemoryCUDA, EstimateMatchesS3X)
{
	const auto dir {temp_dir("gpu_mem_s3x")};
	Background::Background bkg {MultiZone::s3x_background(dir)};
	std::filesystem::remove_all(dir);
	estimate_matches(bkg);
}
