/*
* Tests for Zones::Zones, the description of how the simulation volume is
* split into zones and how they connect.
*/

#include <algorithm>
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
			else if (code == -1) face_code = Zones::face_wall;
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
	// Targets are wall faces and the outer boundary is absorbing, and
	// validate accepts both
	EXPECT_EQ(zones.get_face(2, Zones::ymin), Zones::face_wall);
	EXPECT_EQ(zones.get_face(3, Zones::ymax), Zones::face_wall);
	EXPECT_EQ(zones.get_face(8, Zones::xmax), Zones::face_absorbing);
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

	// Periodic on one side only is allowed, as it always has been for
	// Gkeyll backgrounds
	zones.set_face(0, Zones::ymax, Zones::face_from_bound_type(0, 0));
	zones.validate();
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

	// Same number of cells along the shared face, but different edges
	{
		Zones::Zones zones;
		add_index_zone(zones, 4, 6, 2);
		auto [x, gx] = uniform(3, 0.0, 3.0);
		auto [y, gy] = uniform(6, 0.0, 12.0);
		auto [z, gz] = uniform(2, 0.0, 6.283185307179586);
		zones.add_zone(x, y, z, gx, gy, gz);
		zones.set_face(0, Zones::xmax, 1);
		zones.set_face(1, Zones::xmin, 0);
		EXPECT_EXIT(zones.validate(), testing::ExitedWithCode(EXIT_FAILURE),
			"different y cell edges");
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


// -------------------------------------------------------------------------
// Ghost cells
// -------------------------------------------------------------------------

namespace
{
	// Zones in index space with the given ghost width, the same layout as
	// in test_cross_zones_cpu:
	//
	//        +-------+
	//        |   2   |  ny = 5
	//        +-------+
	//   core |   0   | 1 |  ny = 6
	//        +-------+---+
	//         nx = 4  nx = 3
	//
	// Every zone is periodic in z with nz = 3, and every other face in x
	// and y is absorbing.
	Zones::Zones make_three_zones(int ghost)
	{
		Zones::Zones zones(ghost);
		add_index_zone(zones, 4, 6, 3);
		add_index_zone(zones, 3, 6, 3);
		add_index_zone(zones, 4, 5, 3);
		zones.set_face(0, Zones::xmax, 1);
		zones.set_face(1, Zones::xmin, 0);
		zones.set_face(0, Zones::ymax, 2);
		zones.set_face(2, Zones::ymin, 0);
		zones.set_face(0, Zones::xmin, Zones::face_core);
		for (int z {}; z < 3; ++z)
		{
			zones.set_face(z, Zones::zmin, z);
			zones.set_face(z, Zones::zmax, z);
		}
		zones.validate();
		return zones;
	}

	// Where each zone of make_three_zones sits in a global (x, y) index
	// space
	const std::array<std::array<int, 2>, 3> three_zone_offsets
		{{{0, 0}, {4, 0}, {0, 6}}};

	// A field that is linear in the global indices and periodic in k
	double global_field(int t, int gi, int gj, int k, int nz)
	{
		const int kk {((k % nz) + nz) % nz};
		return 1.0 + 100.0 * t + 2.0 * gi + 3.0 * gj + 5.0 * kk;
	}

	// ZonedVector shaped like zones, with every interior cell set to
	// global_field and every ghost cell set to -1
	Vectors::ZonedVector<double> make_global_field(const Zones::Zones& zones,
		int dim1)
	{
		const int g {zones.get_ghost()};
		Vectors::ZonedVector<double> v(dim1, zones.get_nx_all(),
			zones.get_ny_all(), zones.get_nz_all(), g);
		for (int z {}; z < zones.get_nzones(); ++z)
		for (int t {}; t < dim1; ++t)
		for (int i {-g}; i < zones.get_nx(z) + g; ++i)
		for (int j {-g}; j < zones.get_ny(z) + g; ++j)
		for (int k {-g}; k < zones.get_nz(z) + g; ++k)
		{
			const bool interior {i >= 0 && i < zones.get_nx(z) && j >= 0
				&& j < zones.get_ny(z) && k >= 0 && k < zones.get_nz(z)};
			v(z, t, i, j, k) = interior
				? global_field(t, i + three_zone_offsets[z][0],
					j + three_zone_offsets[z][1], k, zones.get_nz(z))
				: -1.0;
		}
		return v;
	}
}

// Ghost cells past one face shared with another zone (or periodic in z)
// get the value at their position in the global layout, in every time
// frame
TEST(ZonesCPU, FillGhostsEdges)
{
	const Zones::Zones zones {make_three_zones(2)};
	auto v {make_global_field(zones, 2)};
	zones.fill_ghosts(v);

	int checked {};
	for (int z {}; z < 3; ++z)
	for (int t {}; t < 2; ++t)
	for (int i {-2}; i < zones.get_nx(z) + 2; ++i)
	for (int j {-2}; j < zones.get_ny(z) + 2; ++j)
	for (int k {-2}; k < zones.get_nz(z) + 2; ++k)
	{
		const bool out_x {i < 0 || i >= zones.get_nx(z)};
		const bool out_y {j < 0 || j >= zones.get_ny(z)};
		if (out_x && out_y) continue;

		// Only ghosts across a face shared with another zone
		if (out_x && zones.get_face(z, i < 0 ? Zones::xmin : Zones::xmax) < 0)
			continue;
		if (out_y && zones.get_face(z, j < 0 ? Zones::ymin : Zones::ymax) < 0)
			continue;

		EXPECT_DOUBLE_EQ(v(z, t, i, j, k), global_field(t,
			i + three_zone_offsets[z][0], j + three_zone_offsets[z][1], k,
			zones.get_nz(z))) << z << " " << t << " " << i << " " << j
			<< " " << k;
		++checked;
	}
	EXPECT_GT(checked, 0);

	// Interior cells are left alone
	for (int z {}; z < 3; ++z)
	for (int t {}; t < 2; ++t)
	for (int i {}; i < zones.get_nx(z); ++i)
	for (int j {}; j < zones.get_ny(z); ++j)
	for (int k {}; k < zones.get_nz(z); ++k)
	{
		ASSERT_DOUBLE_EQ(v(z, t, i, j, k), global_field(t,
			i + three_zone_offsets[z][0], j + three_zone_offsets[z][1], k,
			zones.get_nz(z)));
	}
}

// Ghost cells past absorbing and core faces copy the zone's edge cell
TEST(ZonesCPU, FillGhostsBoundaryFaces)
{
	const Zones::Zones zones {make_three_zones(2)};
	auto v {make_global_field(zones, 1)};
	zones.fill_ghosts(v);

	// Core face of zone 0
	EXPECT_DOUBLE_EQ(v(0, 0, -1, 2, 1), v(0, 0, 0, 2, 1));
	EXPECT_DOUBLE_EQ(v(0, 0, -2, 2, 1), v(0, 0, 0, 2, 1));

	// Absorbing faces of zone 1
	EXPECT_DOUBLE_EQ(v(1, 0, 4, 3, 0), v(1, 0, 2, 3, 0));
	EXPECT_DOUBLE_EQ(v(1, 0, 1, -1, 0), v(1, 0, 1, 0, 0));
}

// Corner ghosts. Where only one order of crossing reaches another zone,
// that's the source. Where neither does, the zone's own corner cell is.
TEST(ZonesCPU, FillGhostsCorners)
{
	const Zones::Zones zones {make_three_zones(2)};
	auto v {make_global_field(zones, 1)};
	zones.fill_ghosts(v);

	// Zone 1 past xmin and ymax: through zone 0 into zone 2, at global
	// (3, 6). Crossing ymax first hits an absorbing face.
	EXPECT_DOUBLE_EQ(v(1, 0, -1, 6, 1), global_field(0, 3, 6, 1, 3));

	// Zone 2 past xmax and ymin: through zone 0 into zone 1, at global
	// (4, 5)
	EXPECT_DOUBLE_EQ(v(2, 0, 4, -1, 1), global_field(0, 4, 5, 1, 3));

	// Zone 0 past xmax and ymax: both orders hit an absorbing face, so it
	// copies its own corner cell
	EXPECT_DOUBLE_EQ(v(0, 0, 4, 6, 1), v(0, 0, 3, 5, 1));
}

// Where the two orders of crossing reach different cells, as at an
// X-point, the ghost is the average of them:
//
//   C | E
//   --+--
//   A | B   with D above B
//
// Past A's xmax then ymax is D, past A's ymax then xmax is E.
TEST(ZonesCPU, FillGhostsXPointAverage)
{
	Zones::Zones zones(1);
	for (int z {}; z < 5; ++z) add_index_zone(zones, 2, 2, 1);
	const int A {0}, B {1}, C {2}, D {3}, E {4};
	zones.set_face(A, Zones::xmax, B);
	zones.set_face(B, Zones::xmin, A);
	zones.set_face(A, Zones::ymax, C);
	zones.set_face(C, Zones::ymin, A);
	zones.set_face(B, Zones::ymax, D);
	zones.set_face(D, Zones::ymin, B);
	zones.set_face(C, Zones::xmax, E);
	zones.set_face(E, Zones::xmin, C);
	zones.validate();

	EXPECT_EQ(zones.resolve_cell({A, 2, 2, 0}),
		(std::vector<Zones::Cell> {{D, 0, 0, 0}, {E, 0, 0, 0}}));

	Vectors::ZonedVector<double> v(1, zones.get_nx_all(),
		zones.get_ny_all(), zones.get_nz_all(), 1);
	v(D, 0, 0, 0) = 10.0;
	v(E, 0, 0, 0) = 20.0;
	zones.fill_ghosts(v);
	EXPECT_DOUBLE_EQ(v(A, 0, 2, 2, 0), 15.0);
}

// A neighbor thinner than the ghost width is stepped through into the zone
// beyond it
TEST(ZonesCPU, FillGhostsThinNeighbor)
{
	Zones::Zones zones(2);
	add_index_zone(zones, 4, 3, 1);
	add_index_zone(zones, 1, 3, 1);
	add_index_zone(zones, 3, 3, 1);
	zones.set_face(0, Zones::xmax, 1);
	zones.set_face(1, Zones::xmin, 0);
	zones.set_face(1, Zones::xmax, 2);
	zones.set_face(2, Zones::xmin, 1);
	zones.validate();

	EXPECT_EQ(zones.resolve_cell({0, 5, 1, 0}),
		(std::vector<Zones::Cell> {{2, 0, 1, 0}}));
	EXPECT_EQ(zones.resolve_cell({2, -2, 1, 0}),
		(std::vector<Zones::Cell> {{0, 3, 1, 0}}));
}

// Every ghost cell of the real DIII-D mesh resolves to interior cells
TEST(ZonesCPU, FillGhostsDIIIDMesh)
{
	const Zones::Zones zones {make_diiid_zones()};
	Vectors::ZonedVector<double> v(1, zones.get_nx_all(),
		zones.get_ny_all(), zones.get_nz_all(), 2);
	for (int z {}; z < zones.get_nzones(); ++z)
	for (int i {}; i < zones.get_nx(z); ++i)
	for (int j {}; j < zones.get_ny(z); ++j)
	for (int k {}; k < zones.get_nz(z); ++k)
	{
		v(z, i, j, k) = 1.0;
	}
	zones.fill_ghosts(v);
	for (const double val : v.get_data()) ASSERT_DOUBLE_EQ(val, 1.0);
}

// A vector that doesn't match the zones stops the run
TEST(ZonesCPU, FillGhostsMismatch)
{
	const Zones::Zones zones {make_three_zones(2)};
	Vectors::ZonedVector<double> v(1, zones.get_nx_all(),
		zones.get_ny_all(), zones.get_nz_all(), 1);
	EXPECT_EXIT(zones.fill_ghosts(v), testing::ExitedWithCode(EXIT_FAILURE),
		"don't match the zones");
}

// Four zones in a 2x2 block, so every corner where they meet is regular
// and both orders of crossing reach the same cell:
//
//   2 | 3
//   --+--
//   0 | 1
//
// Every ghost cell that lands inside the block, including corner ghosts in
// x and y and in x, y and z, gets the value at its global position.
TEST(ZonesCPU, FillGhostsRegularCorners)
{
	Zones::Zones zones(2);
	const std::array<std::array<int, 2>, 4> offsets {{{0, 0}, {3, 0},
		{0, 3}, {3, 3}}};
	for (int z {}; z < 4; ++z)
	{
		add_index_zone(zones, 3, 3, 3);
		zones.set_face(z, Zones::zmin, z);
		zones.set_face(z, Zones::zmax, z);
	}
	auto connect = [&](int a, Zones::Face face, int b)
	{
		zones.set_face(a, face, b);
		zones.set_face(b, static_cast<Zones::Face>(face ^ 1), a);
	};
	connect(0, Zones::xmax, 1);
	connect(2, Zones::xmax, 3);
	connect(0, Zones::ymax, 2);
	connect(1, Zones::ymax, 3);
	zones.validate();

	EXPECT_EQ(zones.resolve_cell({0, 3, 3, 0}),
		(std::vector<Zones::Cell> {{3, 0, 0, 0}}));
	EXPECT_EQ(zones.resolve_cell({3, -1, -2, 4}),
		(std::vector<Zones::Cell> {{0, 2, 1, 1}}));

	Vectors::ZonedVector<double> v(2, zones.get_nx_all(),
		zones.get_ny_all(), zones.get_nz_all(), 2);
	for (int z {}; z < 4; ++z)
	for (int t {}; t < 2; ++t)
	for (int i {}; i < 3; ++i)
	for (int j {}; j < 3; ++j)
	for (int k {}; k < 3; ++k)
	{
		v(z, t, i, j, k) = global_field(t, i + offsets[z][0],
			j + offsets[z][1], k, 3);
	}
	zones.fill_ghosts(v);

	int corners {};
	for (int z {}; z < 4; ++z)
	for (int t {}; t < 2; ++t)
	for (int i {-2}; i < 5; ++i)
	for (int j {-2}; j < 5; ++j)
	for (int k {-2}; k < 5; ++k)
	{
		const int gi {i + offsets[z][0]};
		const int gj {j + offsets[z][1]};
		if (gi < 0 || gi >= 6 || gj < 0 || gj >= 6) continue;
		EXPECT_DOUBLE_EQ(v(z, t, i, j, k), global_field(t, gi, gj, k, 3))
			<< z << " " << t << " " << i << " " << j << " " << k;
		corners += ((i < 0 || i >= 3) && (j < 0 || j >= 3));
	}
	EXPECT_GT(corners, 0);
}

// Ghost cells wrap around a ring of zones, like the core
TEST(ZonesCPU, FillGhostsRing)
{
	Zones::Zones zones(2);
	add_index_zone(zones, 3, 4, 1);
	add_index_zone(zones, 3, 7, 1);
	zones.set_face(0, Zones::ymax, 1);
	zones.set_face(1, Zones::ymin, 0);
	zones.set_face(1, Zones::ymax, 0);
	zones.set_face(0, Zones::ymin, 1);
	zones.validate();

	Vectors::ZonedVector<double> v(1, zones.get_nx_all(),
		zones.get_ny_all(), zones.get_nz_all(), 2);
	for (int z {}; z < 2; ++z)
	for (int i {}; i < 3; ++i)
	for (int j {}; j < zones.get_ny(z); ++j)
	{
		v(z, i, j, 0) = 100 * z + 10 * i + j;
	}
	zones.fill_ghosts(v);

	for (int i {}; i < 3; ++i)
	{
		EXPECT_DOUBLE_EQ(v(0, i, -1, 0), v(1, i, 6, 0));
		EXPECT_DOUBLE_EQ(v(0, i, -2, 0), v(1, i, 5, 0));
		EXPECT_DOUBLE_EQ(v(0, i, 4, 0), v(1, i, 0, 0));
		EXPECT_DOUBLE_EQ(v(1, i, -1, 0), v(0, i, 3, 0));
		EXPECT_DOUBLE_EQ(v(1, i, 7, 0), v(0, i, 0, 0));
		EXPECT_DOUBLE_EQ(v(1, i, 8, 0), v(0, i, 1, 0));
	}
}

// Changing a face after fill_ghosts has run makes the next call use the
// new connections, not the cached ones
TEST(ZonesCPU, FillGhostsAfterFaceChange)
{
	Zones::Zones zones {make_three_zones(2)};
	auto v {make_global_field(zones, 1)};
	zones.fill_ghosts(v);
	EXPECT_DOUBLE_EQ(v(0, 0, 4, 2, 1), global_field(0, 4, 2, 1, 3));

	zones.set_face(0, Zones::xmax, Zones::face_absorbing);
	zones.set_face(1, Zones::xmin, Zones::face_absorbing);
	zones.fill_ghosts(v);
	EXPECT_DOUBLE_EQ(v(0, 0, 4, 2, 1), v(0, 0, 3, 2, 1));
}

// Spot checks on the real DIII-D mesh (0-based zones): across the core's
// psi+ face into the SOL, around the core ring, and the X-point corner
// averaged from the divertor leg and SOL zones on either side
TEST(ZonesCPU, ResolveDIIIDMesh)
{
	const Zones::Zones zones {make_diiid_zones()};

	EXPECT_EQ(zones.resolve_cell({0, 67, 10, 3}),
		(std::vector<Zones::Cell> {{6, 0, 10, 3}}));
	EXPECT_EQ(zones.resolve_cell({0, 10, 976, 3}),
		(std::vector<Zones::Cell> {{1, 10, 0, 3}}));
	EXPECT_EQ(zones.resolve_cell({1, 10, -1, 3}),
		(std::vector<Zones::Cell> {{0, 10, 975, 3}}));

	auto xpoint {zones.resolve_cell({0, 67, 976, 3})};
	std::ranges::sort(xpoint, {}, &Zones::Cell::zone);
	EXPECT_EQ(xpoint,
		(std::vector<Zones::Cell> {{5, 0, 0, 3}, {7, 0, 0, 3}}));
}
