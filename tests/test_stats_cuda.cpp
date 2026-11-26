/*
* GPU statistics and births on multi-zone backgrounds:
*
* - record_stats_gpu adds each particle to its own zone's cell, the same
*   as record_stats_cpu: the same counts, and sums equal to rounding (the
*   order of the atomic additions differs).
* - fill_slots_gpu starts particles in imp_zone_start_val, inside that
*   zone, with the cell indices of their positions in it.
*/

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "background.h"
#include "background_device.h"
#include "impurity_stats.h"
#include "impurity_stats_device.h"
#include "multizone_cases.h"
#include "options.h"
#include "options_device.h"
#include "slots.h"
#include "slots_device.h"

#include "impurity_stats.cuh"
#include "init_rngs.cuh"
#include "slots.cuh"

namespace
{
	std::filesystem::path test_dir(const std::string& name)
	{
		return std::filesystem::temp_directory_path() / ("flan_test_stats_cuda_"
			+ name + "_" + std::to_string(getpid()));
	}

	void stats_match_cpu(const Background::Background& bkg_in)
	{
		Background::Background bkg {bkg_in};
		Slots::Slots slots {MultiZone::spread_particles(bkg, 300, 11)};

		// Some dead particles, which neither records
		for (int i {}; i < slots.N(); i += 7) slots.set_state(i, 1);
		const Options::Options opts {};
		const double dt {1e-8};

		Impurity::Statistics cpu {bkg};
		Impurity::record_stats_cpu(cpu, slots, opts, dt);

		Impurity::Statistics gpu {bkg};
		ImpurityStats::StatisticsDevice stats_d {gpu.to_device(0)};
		Slots::SlotsDevice slots_d {slots.to_device()};
		ImpurityStats::record_stats_gpu(stats_d, slots_d, dt);
		gpu.add_stats_device(stats_d, 0);
		Slots::free_slots(slots_d);
		Impurity::free_stats(stats_d, 0);

		EXPECT_EQ(gpu.get_counts().get_data(), cpu.get_counts().get_data());
		int recorded {};
		for (const int c : cpu.get_counts().get_data()) recorded += c;
		EXPECT_EQ(recorded, slots.N() - (slots.N() + 6) / 7);

		const std::vector<std::pair<const std::vector<double>*,
			const std::vector<double>*>> sums {
			{&gpu.get_weights().get_data(), &cpu.get_weights().get_data()},
			{&gpu.get_vX().get_data(), &cpu.get_vX().get_data()},
			{&gpu.get_vZ().get_data(), &cpu.get_vZ().get_data()},
			{&gpu.get_vY().get_data(), &cpu.get_vY().get_data()},
			{&gpu.get_charge().get_data(), &cpu.get_charge().get_data()}};
		for (const auto& [g, c] : sums)
		{
			double scale {};
			for (const double v : *c) scale = std::max(scale, std::abs(v));
			ASSERT_GT(scale, 0.0);
			for (std::size_t n {}; n < c->size(); ++n)
				ASSERT_NEAR((*g)[n], (*c)[n], 1e-12 * scale) << "index " << n;
		}
	}
}

TEST(StatsCUDA, MatchesCPUOnS3XGrid)
{
	const auto dir {test_dir("s3x")};
	const Background::Background bkg {MultiZone::s3x_background(dir)};
	std::filesystem::remove_all(dir);
	stats_match_cpu(bkg);
}

TEST(StatsCUDA, MatchesCPUOnHermesGrid)
{
	const auto dir {test_dir("hermes")};
	const Background::Background bkg {MultiZone::hermes_background(dir)};
	std::filesystem::remove_all(dir);
	stats_match_cpu(bkg);
}

// Births over the whole of zone 3 (imp_zone_start_val = 3) of the
// SOLEDGE3X mesh, at the main ion temperature
TEST(StatsCUDA, BirthsInStartZone)
{
	const auto dir {test_dir("births")};
	Background::Background bkg {MultiZone::s3x_background(dir)};
	std::filesystem::remove_all(dir);
	ASSERT_EQ(bkg.get_zones().get_nzones(), 3);

	Options::Options opts {};
	opts.set_imp_zone_start_val(3);
	opts.set_imp_tstart_opt("full_range");
	opts.set_imp_xstart_opt("full_range");
	opts.set_imp_ystart_opt("full_range");
	opts.set_imp_zstart_opt("full_range");
	opts.set_imp_temp_start_opt("main_ion");
	opts.set_imp_num(2000);

	const int n {2000};
	Slots::Slots slots(n);
	for (int i {}; i < n; ++i) slots.set_state(i, 1);
	Slots::SlotsDevice slots_d {slots.to_device()};
	pcg32* rngs_d {init_rngs_cuda(slots_d, 5)};
	Options::OptionsDevice* opts_d {opts.to_device()};
	Background::BackgroundDevice bkg_d {bkg.to_device()};
	int rem {n};
	int alive {};
	Slots::fill_slots_gpu(slots_d, rem, alive, rngs_d, opts_d, bkg_d);
	slots = slots.to_host(slots_d);
	free_rngs(rngs_d, slots_d.device_id);
	Slots::free_slots(slots_d);
	Options::free_opts(opts_d);
	Background::free_bkg(bkg_d);

	EXPECT_EQ(alive, n);
	EXPECT_EQ(rem, 0);
	const int zone {2};
	const auto& gx {bkg.get_grid_x(zone)};
	const auto& gy {bkg.get_grid_y(zone)};
	const auto& gz {bkg.get_grid_z(zone)};
	double min_x {gx.back()};
	double max_x {gx.front()};
	for (int i {}; i < n; ++i)
	{
		ASSERT_EQ(slots.state()[i], 0);
		ASSERT_EQ(slots.zone()[i], zone);
		const double x {slots.x()[i]};
		const double y {slots.y()[i]};
		const double z {slots.z()[i]};
		ASSERT_TRUE(x >= gx.front() && x <= gx.back()) << x;
		ASSERT_TRUE(y >= gy.front() && y <= gy.back()) << y;
		ASSERT_TRUE(z >= gz.front() && z <= gz.back()) << z;
		EXPECT_EQ(slots.xidx()[i], Indices::get_nearest_cell_index_cpu(gx,
			x));
		EXPECT_EQ(slots.yidx()[i], Indices::get_nearest_cell_index_cpu(gy,
			y));
		EXPECT_EQ(slots.zidx()[i], Indices::get_nearest_cell_index_cpu(gz,
			z));
		EXPECT_TRUE(std::isfinite(slots.vX()[i]) && slots.vX()[i] != 0.0);
		min_x = std::min(min_x, x);
		max_x = std::max(max_x, x);
	}

	// Spread over the zone, not one corner of it
	const double span {gx.back() - gx.front()};
	EXPECT_LT(min_x, gx.front() + 0.1 * span);
	EXPECT_GT(max_x, gx.back() - 0.1 * span);
}
