/*
* Tests for Impurity::Statistics with more than one zone: recording into
* the right zone and cell, the density from each zone's cell volumes, the
* averages, and packing for MPI.
*/

#include <gtest/gtest.h>
#include <vector>

#include "background.h"
#include "impurity_stats.h"
#include "options.h"
#include "slots.h"
#include "vectors.h"
#include "zones.h"


namespace
{
	// Grid with n cells between lo and hi, returned as (centers, edges)
	std::pair<std::vector<double>, std::vector<double>> uniform(int n,
		double lo, double hi)
	{
		std::vector<double> centers(n), edges(n + 1);
		for (int i {}; i <= n; ++i) edges[i] = lo + i * (hi - lo) / n;
		for (int i {}; i < n; ++i)
			centers[i] = 0.5 * (edges[i] + edges[i + 1]);
		return {centers, edges};
	}

	// Background with two zones of different cell sizes, no ghost cells,
	// two time frames and J = 2 everywhere:
	// - zone 0: 3 x 4 x 2 cells of size 1 x 0.5 x 0.25
	// - zone 1: 2 x 4 x 2 cells of size 2 x 0.5 x 0.25
	Background::Background make_background()
	{
		Zones::Zones zones {};
		auto [x0, gx0] = uniform(3, 0.0, 3.0);
		auto [x1, gx1] = uniform(2, 0.0, 4.0);
		auto [y, gy] = uniform(4, 0.0, 2.0);
		auto [z, gz] = uniform(2, 0.0, 0.5);
		zones.add_zone(x0, y, z, gx0, gy, gz);
		zones.add_zone(x1, y, z, gx1, gy, gz);
		zones.set_face(0, Zones::xmax, 1);
		zones.set_face(1, Zones::xmin, 0);

		Background::Background bkg {};
		std::vector<double> times {0.0, 1.0};
		bkg.move_into_times(times);
		bkg.move_into_zones(zones);

		Vectors::ZonedVector<BkgFPType> J(1, {3, 2}, {4, 4}, {2, 2});
		for (auto& val : J.get_data()) val = 2.0;
		bkg.move_into_J(J);

		// Any time-dependent array sets the number of times
		Vectors::ZonedVector<BkgFPType> ne(2, {3, 2}, {4, 4}, {2, 2});
		bkg.move_into_ne(ne);
		return bkg;
	}
}

// Statistics made from a background have its zones and number of times
TEST(StatsCPU, ShapeFromBackground)
{
	const auto bkg {make_background()};
	Impurity::Statistics stats {bkg};
	EXPECT_EQ(stats.get_dim1(), 2);
	EXPECT_EQ(stats.get_nzones(), 2);
	EXPECT_EQ(stats.get_counts().get_nx(1), 2);
	EXPECT_EQ(stats.get_weights().get_ny(0), 4);
	EXPECT_EQ(stats.get_counts().get_ghost(), 0);
	EXPECT_EQ(stats.get_counts().get_data().size(), 2u * (24 + 16));

	// dim2 to dim4 only exist for a single zone
	EXPECT_EXIT(stats.get_dim2(), testing::ExitedWithCode(EXIT_FAILURE),
		"single zone");
	EXPECT_TRUE(stats.same_shape().get_s().check_same_shape(stats.get_s()));
}

// A particle is recorded in its zone, time frame and cell, and nowhere else
TEST(StatsCPU, RecordInZone)
{
	const auto bkg {make_background()};
	Impurity::Statistics stats {bkg};

	Slots::Slots slots(2);
	for (int i {}; i < 2; ++i)
	{
		slots.set_zone(i, 1);
		slots.set_tidx(i, 1);
		slots.set_xidx(i, 0);
		slots.set_yidx(i, 2);
		slots.set_zidx(i, 1);
		slots.set_vX(i, 100.0);
		slots.set_q(i, 3);
	}
	slots.set_state(1, 1);	// Dead, not recorded

	const Options::Options opts {};
	const double dt {1e-8};
	Impurity::record_stats_cpu(stats, slots, opts, dt);

	const double w_dt {slots.weight()[0] * dt};
	EXPECT_EQ(stats.get_counts()(1, 1, 0, 2, 1), 1);
	EXPECT_DOUBLE_EQ(stats.get_weights()(1, 1, 0, 2, 1), w_dt);
	EXPECT_DOUBLE_EQ(stats.get_vX()(1, 1, 0, 2, 1), 100.0 * w_dt);
	EXPECT_DOUBLE_EQ(stats.get_charge()(1, 1, 0, 2, 1), 3.0 * w_dt);

	int total {};
	for (const int c : stats.get_counts().get_data()) total += c;
	EXPECT_EQ(total, 1);
}

// The density is the weight over each zone's cell volume J dx dy dz
TEST(StatsCPU, DensityPerZone)
{
	const auto bkg {make_background()};
	Impurity::Statistics stats {bkg};
	stats.get_weights()(0, 1, 2, 3, 1) = 6.0;
	stats.get_weights()(1, 0, 1, 0, 0) = 6.0;
	stats.calc_density(bkg, 1);

	EXPECT_DOUBLE_EQ(stats.get_density()(0, 1, 2, 3, 1),
		6.0 / (2.0 * 1.0 * 0.5 * 0.25));
	EXPECT_DOUBLE_EQ(stats.get_density()(1, 0, 1, 0, 0),
		6.0 / (2.0 * 2.0 * 0.5 * 0.25));
	EXPECT_DOUBLE_EQ(stats.get_density()(1, 1, 1, 0, 0), 0.0);
}

// The averages divide the running sums by the weight in cells with counts
// and are zero elsewhere, in every zone
TEST(StatsCPU, Averages)
{
	const auto bkg {make_background()};
	Impurity::Statistics stats {bkg};
	for (int zone {}; zone < 2; ++zone)
	{
		stats.get_counts()(zone, 0, 1, 1, 1) = 2;
		stats.get_weights()(zone, 0, 1, 1, 1) = 4.0;
		stats.get_vY()(zone, 0, 1, 1, 1) = 8.0;
		stats.get_vz()(zone, 0, 1, 1, 1) = -2.0;
		stats.get_charge()(zone, 0, 1, 1, 1) = 12.0;
		stats.get_s()(zone, 0, 1, 1, 1) = 1.0;

		// No counts, so set to zero
		stats.get_vX()(zone, 1, 0, 0, 0) = 5.0;
		stats.get_charge()(zone, 1, 0, 0, 0) = 5.0;
	}
	stats.calc_vels();
	stats.calc_charge();
	stats.calc_s();

	for (int zone {}; zone < 2; ++zone)
	{
		EXPECT_DOUBLE_EQ(stats.get_vY()(zone, 0, 1, 1, 1), 2.0);
		EXPECT_DOUBLE_EQ(stats.get_vz()(zone, 0, 1, 1, 1), -0.5);
		EXPECT_DOUBLE_EQ(stats.get_charge()(zone, 0, 1, 1, 1), 3.0);
		EXPECT_DOUBLE_EQ(stats.get_s()(zone, 0, 1, 1, 1), 0.25);
		EXPECT_DOUBLE_EQ(stats.get_vX()(zone, 1, 0, 0, 0), 0.0);
		EXPECT_DOUBLE_EQ(stats.get_charge()(zone, 1, 0, 0, 0), 0.0);
	}
}

// Packing and unpacking keeps every value of every zone
TEST(StatsCPU, PackUnpack)
{
	Impurity::Statistics stats {2, {3, 2}, {4, 5}, {2, 1}};
	int n {};
	for (auto& c : stats.get_counts().get_data()) c = n++;
	for (auto& w : stats.get_weights().get_data()) w = 0.5 * n++;
	for (auto& v : stats.get_vz().get_data()) v = -1.0 * n++;

	std::vector<double> buf {};
	stats.pack(buf);
	EXPECT_EQ(buf.size(), 10 * stats.get_counts().get_data().size());

	Impurity::Statistics copy {stats.same_shape()};
	copy.unpack(buf);
	EXPECT_EQ(copy.get_counts().get_data(), stats.get_counts().get_data());
	EXPECT_EQ(copy.get_weights().get_data(), stats.get_weights().get_data());
	EXPECT_EQ(copy.get_vz().get_data(), stats.get_vz().get_data());
}
