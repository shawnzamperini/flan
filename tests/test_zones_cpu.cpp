/*
* Tests for Zones::Zones, the description of how the simulation volume is
* split into zones and how they connect.
*/

#include <array>
#include <gtest/gtest.h>
#include <vector>

#include "slots.h"
#include "zones.h"


// Uniform grid with n cells between lo and hi, returned as (centers, edges)
static std::pair<std::vector<double>, std::vector<double>> uniform(int n,
	double lo, double hi)
{
	std::vector<double> centers(n), edges(n + 1);
	const double d {(hi - lo) / n};
	for (int i {}; i <= n; ++i) edges[i] = lo + i * d;
	for (int i {}; i < n; ++i) centers[i] = lo + (i + 0.5) * d;
	return {centers, edges};
}

// Adds a zone that uses index space in x and y, i.e., the cell edges are
// 0, 1, ..., n, the same as the SOLEDGE3X zones will.
static int add_index_zone(Zones::Zones& zones, int nx, int ny, int nz)
{
	auto [x, gx] = uniform(nx, 0.0, nx);
	auto [y, gy] = uniform(ny, 0.0, ny);
	auto [z, gz] = uniform(nz, 0.0, 6.283185307179586);
	return zones.add_zone(x, y, z, gx, gy, gz);
}

// The zones of the DIII-D 167196 SOLEDGE3X mesh (mesh_S3X.h5), with the
// real cell counts and neighbors table. The neighbors codes are 1-based zone
// numbers, -1 = target, -2 = core and -3 = outer boundary.
static Zones::Zones make_diiid_zones()
{
	const std::array<int, 13> npsi {67, 67, 13, 13, 34, 34, 34, 34, 4, 4,
		31, 31, 31};
	const std::array<int, 13> ntheta {976, 648, 98, 53, 98, 53, 976, 648, 98,
		648, 976, 53, 2};
	const std::array<std::array<int, 13>, 4> neighbors {{
		{-2, -2, -3, -3, 3, 4, 1, 2, 5, 8, 7, 6, -3},
		{7, 8, 5, 6, 9, 12, 11, 10, -3, -3, -3, -3, -3},
		{2, 1, -1, 3, -1, 7, 8, 5, -1, 9, 13, 11, -1},
		{2, 1, 4, -1, 8, -1, 6, 7, 10, -1, 12, -1, 11}}};

	Zones::Zones zones(2);
	for (int z {}; z < 13; ++z) add_index_zone(zones, npsi[z], ntheta[z], 16);
	for (int z {}; z < 13; ++z)
	{
		for (int f {}; f < 4; ++f)
		{
			const int code {neighbors[f][z]};
			int face_code {Zones::face_absorbing};
			if (code > 0) face_code = code - 1;
			else if (code == -2) face_code = Zones::face_core;
			zones.set_face(z, static_cast<Zones::Face>(f), face_code);
		}

		// Toroidally periodic
		zones.set_face(z, Zones::zmin, z);
		zones.set_face(z, Zones::zmax, z);
	}
	return zones;
}

// The real SOLEDGE3X connectivity is consistent, and the core ring is zones
// 0 and 1
TEST(ZonesCPU, DIIIDMeshValidates)
{
	Zones::Zones zones {make_diiid_zones()};
	zones.validate();

	EXPECT_EQ(zones.get_nzones(), 13);
	EXPECT_EQ(zones.get_ghost(), 2);
	EXPECT_EQ(zones.get_core_zones(), (std::vector<int> {0, 1}));

	// Core ring: the two zones wrap around into each other in theta
	EXPECT_EQ(zones.get_face(0, Zones::ymin), 1);
	EXPECT_EQ(zones.get_face(0, Zones::ymax), 1);
	EXPECT_EQ(zones.get_face(0, Zones::xmax), 6);
	EXPECT_EQ(zones.get_face(6, Zones::xmin), 0);
	EXPECT_EQ(zones.get_face(2, Zones::ymin), Zones::face_absorbing);
}

TEST(ZonesCPU, Accessors)
{
	Zones::Zones zones {make_diiid_zones()};

	EXPECT_EQ(zones.get_nx(0), 67);
	EXPECT_EQ(zones.get_ny(1), 648);
	EXPECT_EQ(zones.get_nz(12), 16);
	EXPECT_EQ(zones.get_nx_all().size(), 13u);
	EXPECT_EQ(zones.get_ny_all()[3], 53);
	EXPECT_EQ(zones.get_nz_all()[5], 16);
	EXPECT_DOUBLE_EQ(zones.get_x_min(4), 0.0);
	EXPECT_DOUBLE_EQ(zones.get_x_max(4), 34.0);
	EXPECT_DOUBLE_EQ(zones.get_y_max(4), 98.0);
	EXPECT_DOUBLE_EQ(zones.get_z_max(4), 6.283185307179586);
	EXPECT_EQ(zones.get_grid_x(4).size(), 35u);
	EXPECT_DOUBLE_EQ(zones.get_x(4)[0], 0.5);
}

// A Gkeyll-style single zone with faces from the *bound_type options
TEST(ZonesCPU, SingleZoneFromBoundTypes)
{
	EXPECT_EQ(Zones::face_from_bound_type(0, 0), Zones::face_absorbing);
	EXPECT_EQ(Zones::face_from_bound_type(1, 0), 0);
	EXPECT_EQ(Zones::face_from_bound_type(2, 0), Zones::face_core);
	EXPECT_EXIT(Zones::face_from_bound_type(3, 0),
		testing::ExitedWithCode(EXIT_FAILURE), "Bound type 3");

	auto [x, gx] = uniform(10, -0.1, 0.1);
	auto [y, gy] = uniform(20, 0.0, 1.0);
	auto [z, gz] = uniform(5, -1.0, 1.0);
	Zones::Zones zones;
	EXPECT_EQ(zones.add_zone(x, y, z, gx, gy, gz), 0);
	zones.set_face(0, Zones::xmin, Zones::face_from_bound_type(2, 0));
	zones.set_face(0, Zones::ymin, Zones::face_from_bound_type(1, 0));
	zones.set_face(0, Zones::ymax, Zones::face_from_bound_type(1, 0));
	zones.validate();

	EXPECT_EQ(zones.get_ghost(), 0);
	EXPECT_EQ(zones.get_face(0, Zones::xmax), Zones::face_absorbing);
	EXPECT_EQ(zones.get_core_zones(), (std::vector<int> {0}));
	EXPECT_DOUBLE_EQ(zones.get_z_min(0), -1.0);
}

// Each kind of inconsistency stops the run
TEST(ZonesCPU, ValidateFailures)
{
	// Only one side of a shared face points at the other zone
	{
		Zones::Zones zones;
		add_index_zone(zones, 4, 6, 2);
		add_index_zone(zones, 3, 6, 2);
		zones.set_face(0, Zones::xmax, 1);
		EXPECT_EXIT(zones.validate(), testing::ExitedWithCode(EXIT_FAILURE),
			"connects to zone 1");
	}

	// Different number of cells along the shared face
	{
		Zones::Zones zones;
		add_index_zone(zones, 4, 6, 2);
		add_index_zone(zones, 3, 5, 2);
		zones.set_face(0, Zones::xmax, 1);
		zones.set_face(1, Zones::xmin, 0);
		EXPECT_EXIT(zones.validate(), testing::ExitedWithCode(EXIT_FAILURE),
			"different numbers of cells in y");
	}

	// Zone that doesn't exist
	{
		Zones::Zones zones;
		add_index_zone(zones, 4, 6, 2);
		zones.set_face(0, Zones::ymax, 3);
		EXPECT_EXIT(zones.validate(), testing::ExitedWithCode(EXIT_FAILURE),
			"invalid code 3");
	}

	// Edges must be increasing and one longer than the centers
	{
		Zones::Zones zones;
		zones.add_zone({0.5, 1.5}, {0.5}, {0.5}, {0.0, 2.0, 1.0},
			{0.0, 1.0}, {0.0, 1.0});
		EXPECT_EXIT(zones.validate(), testing::ExitedWithCode(EXIT_FAILURE),
			"not increasing");
	}
	{
		Zones::Zones zones;
		zones.add_zone({0.5, 1.5}, {0.5}, {0.5}, {0.0, 1.0}, {0.0, 1.0},
			{0.0, 1.0});
		EXPECT_EXIT(zones.validate(), testing::ExitedWithCode(EXIT_FAILURE),
			"2 x cell centers and 2 edges");
	}
}

// New particle slots start in zone 0, and the zone can be set per slot
TEST(ZonesCPU, SlotsZone)
{
	Slots::Slots slots(4);
	EXPECT_EQ(slots.zone(), (std::vector<int> {0, 0, 0, 0}));
	slots.set_zone(2, 7);
	EXPECT_EQ(slots.zone()[2], 7);
	EXPECT_EQ(slots.zone()[1], 0);
}
