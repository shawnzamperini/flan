/*
* Tests for interpolation across zones: ghost cell centers, the neighbor and
* bracketing cells near zone edges, and fields interpolated by
* Background::interp_ne and Boris::update_velocity_cpu staying continuous
* across a face shared by two zones.
*/

#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <utility>
#include <vector>

#include "background.h"
#include "boris.h"
#include "constants.h"
#include "indices.h"
#include "options.h"
#include "pcg32.h"
#include "slots.h"
#include "utilities.h"
#include "vectors.h"
#include "zones.h"


namespace
{
	constexpr double two_pi {6.283185307179586};

	// Grid from its edges, returned as (centers, edges)
	std::pair<std::vector<double>, std::vector<double>> from_edges(
		std::vector<double> edges)
	{
		std::vector<double> centers(edges.size() - 1);
		for (std::size_t i {}; i < centers.size(); ++i)
			centers[i] = 0.5 * (edges[i] + edges[i + 1]);
		return {centers, edges};
	}

	// Grid with n cells between lo and hi, returned as (centers, edges)
	std::pair<std::vector<double>, std::vector<double>> uniform(int n,
		double lo, double hi)
	{
		std::vector<double> edges(n + 1);
		for (int i {}; i <= n; ++i) edges[i] = lo + i * (hi - lo) / n;
		return from_edges(edges);
	}

	// Zone in index space in x and y (edges 0, 1, ..., n), like SOLEDGE3X,
	// and periodic in z over [0, 2pi]
	int add_index_zone(Zones::Zones& zones, int nx, int ny, int nz)
	{
		auto [x, gx] = uniform(nx, 0.0, nx);
		auto [y, gy] = uniform(ny, 0.0, ny);
		auto [z, gz] = uniform(nz, 0.0, two_pi);
		const int zone {zones.add_zone(x, y, z, gx, gy, gz)};
		zones.set_face(zone, Zones::zmin, zone);
		zones.set_face(zone, Zones::zmax, zone);
		return zone;
	}

	// Two zones side by side in x with the given ghost width. Zone 0 has
	// nx = 4 and zone 1 nx = 3. Both have ny = 3 and nz = 2.
	//
	//   | 0 | 1 |
	//
	// The x and y faces that aren't shared are absorbing.
	Zones::Zones make_two_zones(int ghost)
	{
		Zones::Zones zones(ghost);
		add_index_zone(zones, 4, 3, 2);
		add_index_zone(zones, 3, 3, 2);
		zones.set_face(0, Zones::xmax, 1);
		zones.set_face(1, Zones::xmin, 0);
		zones.validate();
		return zones;
	}

	// Offset of each zone of make_two_zones in a global x index space
	constexpr std::array<double, 2> x_offsets {0.0, 4.0};

	// Field that is linear in the global x and in y, and the same in every
	// time frame and toroidal plane
	double linear(double gx, double y) {return 1.0 + 2.0 * gx + 3.0 * y;}

	// ZonedVector shaped like zones with dim1 time frames, set to f(zone,
	// t, x, y) at the interior cell centers, with the ghost cells then
	// filled
	template <typename F>
	Vectors::ZonedVector<BkgFPType> make_field_t(const Zones::Zones& zones,
		int dim1, F f)
	{
		Vectors::ZonedVector<BkgFPType> v(dim1, zones.get_nx_all(),
			zones.get_ny_all(), zones.get_nz_all(), zones.get_ghost());
		for (int z {}; z < zones.get_nzones(); ++z)
		for (int t {}; t < dim1; ++t)
		for (int i {}; i < zones.get_nx(z); ++i)
		for (int j {}; j < zones.get_ny(z); ++j)
		for (int k {}; k < zones.get_nz(z); ++k)
		{
			v(z, t, i, j, k) = f(z, t, zones.get_x(z)[i],
				zones.get_y(z)[j]);
		}
		zones.fill_ghosts(v);
		return v;
	}

	// The same, for an f(zone, x, y) that doesn't depend on time
	template <typename F>
	Vectors::ZonedVector<BkgFPType> make_field(const Zones::Zones& zones,
		int dim1, F f)
	{
		return make_field_t(zones, dim1, [&](int z, int, double x, double y)
			{return f(z, x, y);});
	}

	// Background over make_two_zones(2) with times {0, 1}
	Background::Background make_background()
	{
		Background::Background bkg {};
		std::vector<double> times {0.0, 1.0};
		bkg.move_into_times(times);
		auto zones {make_two_zones(2)};
		bkg.move_into_zones(zones);
		return bkg;
	}
}

// Ghost cell centers in evenly spaced zones are the next centers along
TEST(ZoneInterpCPU, CenterIndexSpace)
{
	const auto zones {make_two_zones(2)};
	EXPECT_DOUBLE_EQ(zones.center(0, 0, 2), 2.5);
	EXPECT_DOUBLE_EQ(zones.center(0, 0, 4), 4.5);
	EXPECT_DOUBLE_EQ(zones.center(0, 0, 5), 5.5);
	EXPECT_DOUBLE_EQ(zones.center(1, 0, -1), -0.5);
	EXPECT_DOUBLE_EQ(zones.center(1, 0, -2), -1.5);

	// Periodic in z
	EXPECT_DOUBLE_EQ(zones.center(0, 2, -1), 0.75 * two_pi - two_pi);
	EXPECT_DOUBLE_EQ(zones.center(0, 2, 2), 0.25 * two_pi + two_pi);
}

// With different cell widths, a ghost center is the neighbor's center
// shifted into this zone's frame
TEST(ZoneInterpCPU, CenterUnevenCells)
{
	Zones::Zones zones(1);
	auto [x0, gx0] = from_edges({0.0, 1.0, 3.0});
	auto [x1, gx1] = from_edges({0.0, 0.5, 1.0});
	auto [y, gy] = uniform(2, 0.0, 2.0);
	auto [z, gz] = uniform(1, 0.0, two_pi);
	zones.add_zone(x0, y, z, gx0, gy, gz);
	zones.add_zone(x1, y, z, gx1, gy, gz);
	zones.set_face(0, Zones::xmax, 1);
	zones.set_face(1, Zones::xmin, 0);
	zones.validate();

	EXPECT_DOUBLE_EQ(zones.center(0, 0, 2), 3.25);
	EXPECT_DOUBLE_EQ(zones.center(1, 0, -1), -1.0);
}

// A ghost center beyond a neighbor thinner than the ghost width
TEST(ZoneInterpCPU, CenterThinNeighbor)
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

	EXPECT_DOUBLE_EQ(zones.center(0, 0, 5), 5.5);
	EXPECT_DOUBLE_EQ(zones.center(2, 0, -2), -1.5);
}

// There's no ghost center past an absorbing face
TEST(ZoneInterpCPU, CenterPastAbsorbing)
{
	const auto zones {make_two_zones(2)};
	EXPECT_EXIT(zones.center(1, 0, 3), testing::ExitedWithCode(EXIT_FAILURE),
		"no cell center");
}

// Without ghost cells, or away from shared faces, neighbor_index and
// bracket are exactly get_neighbor_index and bracket_indices
TEST(ZoneInterpCPU, NeighborAndBracketMatchUtilities)
{
	pcg32 rng(11, 0);
	for (const int ghost : {0, 2})
	{
		const auto zones {make_two_zones(ghost)};
		for (int n {}; n < 2000; ++n)
		{
			const int zone {static_cast<int>(rng.next_double() * 2)};
			const int dim {static_cast<int>(rng.next_double() * 3)};
			const double lo {zones.get_min(zone, dim)};
			const double hi {zones.get_max(zone, dim)};
			const double coord {lo + rng.next_double() * (hi - lo)};
			const auto& centers {zones.get_centers(zone, dim)};
			const auto& grid {dim == 0 ? zones.get_grid_x(zone)
				: dim == 1 ? zones.get_grid_y(zone) : zones.get_grid_z(zone)};
			const int idx {Indices::get_nearest_cell_index_cpu(grid, coord)};

			// The faces that could use ghost cells are zone 0's xmax, zone
			// 1's xmin and z, so skip those edges when there are ghosts
			const bool edge_cell {idx == 0
				|| idx == static_cast<int>(centers.size()) - 1};
			const bool shared {dim == 2 || (dim == 0 && edge_cell)};
			if (ghost > 0 && shared) continue;

			EXPECT_EQ(zones.neighbor_index(zone, dim, coord, idx),
				Utilities::get_neighbor_index(coord, centers, idx));
			EXPECT_EQ(zones.bracket(zone, dim, coord),
				Utilities::bracket_indices(centers, coord));
		}
	}
}

// At a shared face with ghost cells, the neighbor and bracketing cells
// include the ghost cell. At an absorbing face they don't.
TEST(ZoneInterpCPU, NeighborAndBracketUseGhosts)
{
	const auto zones {make_two_zones(2)};

	// Zone 0, last cell (center 3.5), outer side: ghost 4
	EXPECT_EQ(zones.neighbor_index(0, 0, 3.8, 3), 4);
	EXPECT_EQ(zones.bracket(0, 0, 3.8), std::make_pair(3, 4));

	// Inner side of the same cell: interior neighbor
	EXPECT_EQ(zones.neighbor_index(0, 0, 3.2, 3), 2);

	// Zone 1, first cell (center 0.5), outer side: ghost -1
	EXPECT_EQ(zones.neighbor_index(1, 0, 0.1, 0), -1);
	EXPECT_EQ(zones.bracket(1, 0, 0.1), std::make_pair(-1, 0));

	// Zone 0's xmin is absorbing, so it holds at the edge cells
	EXPECT_EQ(zones.neighbor_index(0, 0, 0.2, 0), 1);
	EXPECT_EQ(zones.bracket(0, 0, 0.2), std::make_pair(0, 1));

	// Periodic z wraps through the ghost cells
	EXPECT_EQ(zones.neighbor_index(0, 2, 0.1, 0), -1);
	EXPECT_EQ(zones.bracket(0, 2, two_pi - 0.1), std::make_pair(1, 2));
}

// interp_ne is continuous across the face shared by the two zones. Points
// near the face in either zone get the exact value of the linear field at
// their global position. Without the ghost cells, those past the last cell
// center would be held at the edge cell's value.
TEST(ZoneInterpCPU, InterpContinuousAcrossZones)
{
	auto bkg {make_background()};
	auto ne {make_field(bkg.get_zones(), 2, [](int z, double x, double y)
		{return linear(x + x_offsets[z], y);})};
	bkg.move_into_ne(ne);

	const double y {1.2};
	for (const double x : {3.2, 3.6, 3.9, 4.0})
	{
		EXPECT_NEAR(bkg.interp_ne(0, 0.5, x, y, 1.0), linear(x, y), 1e-12)
			<< "zone 0, x = " << x;
	}
	for (const double x : {0.0, 0.1, 0.4, 0.8})
	{
		EXPECT_NEAR(bkg.interp_ne(1, 0.5, x, y, 1.0), linear(x + 4.0, y),
			1e-12) << "zone 1, x = " << x;
	}

	// The absorbing face at zone 0's xmin still holds at the edge cell
	EXPECT_NEAR(bkg.interp_ne(0, 0.5, 0.2, y, 1.0), linear(0.5, y), 1e-12);
}

// Boris sees the same continuous field. With B = 0 and E = (Ex, 0, 0)
// linear in x and time, a particle's velocity changes by q E / m dt using
// the exact local Ex, on both sides of the shared face.
TEST(ZoneInterpCPU, BorisContinuousAcrossZones)
{
	auto bkg {make_background()};
	const auto& zones {bkg.get_zones()};
	// E is 10 V/m larger in the second time frame, so the particle at
	// t = 0.5 sees 5 V/m more than in the first
	auto global_ex = [](int z, int t, double x, double y)
		{return linear(x + x_offsets[z], y) + 10.0 * t;};
	auto global_lin = [](int z, double x, double y)
		{return linear(x + x_offsets[z], y);};
	auto zero = [](int, double, double) {return 0.0;};
	auto one = [](int, double, double) {return 1.0;};

	auto eX {make_field_t(zones, 2, global_ex)};
	auto emag {make_field_t(zones, 2, global_ex)};
	auto eY {make_field(zones, 2, zero)};
	auto eZ {make_field(zones, 2, zero)};
	auto bX {make_field(zones, 2, zero)};
	auto bY {make_field(zones, 2, zero)};
	auto bZ {make_field(zones, 2, zero)};
	auto bmag {make_field(zones, 2, zero)};
	bkg.move_into_eX(eX);
	bkg.move_into_eY(eY);
	bkg.move_into_eZ(eZ);
	bkg.move_into_emag(emag);
	bkg.move_into_bX(bX);
	bkg.move_into_bY(bY);
	bkg.move_into_bZ(bZ);
	bkg.move_into_bmag(bmag);

	// Reciprocal basis that is the identity except for dxdX, which is
	// linear in x like E, so the computational velocity vx = dxdX vX also
	// checks the basis is interpolated through the ghost cells
	auto dxdX {make_field(zones, 1, global_lin)};
	auto dydY {make_field(zones, 1, one)};
	auto dzdZ {make_field(zones, 1, one)};
	auto dxdY {make_field(zones, 1, zero)};
	auto dxdZ {make_field(zones, 1, zero)};
	auto dydX {make_field(zones, 1, zero)};
	auto dydZ {make_field(zones, 1, zero)};
	auto dzdX {make_field(zones, 1, zero)};
	auto dzdY {make_field(zones, 1, zero)};
	bkg.move_into_dxdX(dxdX);
	bkg.move_into_dxdY(dxdY);
	bkg.move_into_dxdZ(dxdZ);
	bkg.move_into_dydX(dydX);
	bkg.move_into_dydY(dydY);
	bkg.move_into_dydZ(dydZ);
	bkg.move_into_dzdX(dzdX);
	bkg.move_into_dzdY(dzdY);
	bkg.move_into_dzdZ(dzdZ);

	// Particles near the face in each zone
	const std::vector<std::pair<int, double>> particles {{0, 3.6}, {0, 3.9},
		{1, 0.1}, {1, 0.4}};
	const double y {1.2};
	const double z {1.0};
	const int n {static_cast<int>(particles.size())};
	Slots::Slots slots(n);
	slots.set_mass(12.0 * Constants::amu_to_kg);
	for (int i {}; i < n; ++i)
	{
		const auto [zone, x] = particles[i];
		slots.set_zone(i, zone);
		slots.set_t(i, 0.5);
		slots.set_x(i, x);
		slots.set_y(i, y);
		slots.set_z(i, z);
		slots.set_q(i, 1);
		slots.set_tidx(i, 0);
		slots.set_xidx(i, Indices::get_nearest_cell_index_cpu(
			zones.get_grid_x(zone), x));
		slots.set_yidx(i, Indices::get_nearest_cell_index_cpu(
			zones.get_grid_y(zone), y));
		slots.set_zidx(i, Indices::get_nearest_cell_index_cpu(
			zones.get_grid_z(zone), z));
	}

	const Options::Options opts {};
	const double dt {1e-9};
	Boris::update_velocity_cpu(slots, bkg, opts, dt);

	const double q_m {1 * -Constants::charge_e / slots.mass()};
	for (int i {}; i < n; ++i)
	{
		const auto [zone, x] = particles[i];
		const double gx {x + x_offsets[zone]};
		const double want {q_m * (linear(gx, y) + 5.0) * dt};
		EXPECT_NEAR(slots.vX()[i], want, 1e-12 * std::abs(want))
			<< "zone " << zone << ", x = " << x;
		EXPECT_NEAR(slots.vY()[i], 0.0, 1e-12 * std::abs(want));
		EXPECT_NEAR(slots.vx()[i], linear(gx, y) * slots.vX()[i],
			1e-12 * std::abs(linear(gx, y) * want))
			<< "zone " << zone << ", x = " << x;
	}
}

// Data loaded as a ZonedVector must match the zones
TEST(ZoneInterpCPU, ZonedSettersCheckShape)
{
	auto bkg {make_background()};
	const auto& zones {bkg.get_zones()};

	// Wrong ghost width
	Vectors::ZonedVector<BkgFPType> no_ghosts(2, zones.get_nx_all(),
		zones.get_ny_all(), zones.get_nz_all(), 0);
	EXPECT_EXIT(bkg.move_into_ne(no_ghosts),
		testing::ExitedWithCode(EXIT_FAILURE), "ghost width");

	// Wrong number of cells
	Vectors::ZonedVector<BkgFPType> wrong_cells(2, {4, 4},
		zones.get_ny_all(), zones.get_nz_all(), 2);
	EXPECT_EXIT(bkg.move_into_te(wrong_cells),
		testing::ExitedWithCode(EXIT_FAILURE), "wrong number of cells");

	// Geometry must have one entry in the leading dimension
	Vectors::ZonedVector<BkgFPType> two_frames(2, zones.get_nx_all(),
		zones.get_ny_all(), zones.get_nz_all(), 2);
	EXPECT_EXIT(bkg.move_into_J(two_frames),
		testing::ExitedWithCode(EXIT_FAILURE), "should be 1");

	// Different number of times than data already loaded
	Vectors::ZonedVector<BkgFPType> ne(2, zones.get_nx_all(),
		zones.get_ny_all(), zones.get_nz_all(), 2);
	bkg.move_into_ne(ne);
	Vectors::ZonedVector<BkgFPType> three_frames(3, zones.get_nx_all(),
		zones.get_ny_all(), zones.get_nz_all(), 2);
	EXPECT_EXIT(bkg.move_into_ti(three_frames),
		testing::ExitedWithCode(EXIT_FAILURE), "dim1");
}

// A field that changes between time frames is interpolated in time on top
// of across the zones
TEST(ZoneInterpCPU, InterpTimeAcrossZones)
{
	auto bkg {make_background()};
	auto ne {make_field_t(bkg.get_zones(), 2,
		[](int z, int t, double x, double y)
		{return linear(x + x_offsets[z], y) + 10.0 * t;})};
	bkg.move_into_ne(ne);

	EXPECT_NEAR(bkg.interp_ne(0, 0.25, 3.9, 1.2, 1.0),
		linear(3.9, 1.2) + 2.5, 1e-12);
	EXPECT_NEAR(bkg.interp_ne(1, 0.75, 0.1, 1.2, 1.0),
		linear(4.1, 1.2) + 7.5, 1e-12);
}

// Near the corner where four zones meet in a 2x2 block, interpolation
// uses the corner ghosts and is still exact for a field linear in x and y,
// from every zone
//
//   2 | 3
//   --+--
//   0 | 1
TEST(ZoneInterpCPU, InterpAcrossCorner)
{
	Zones::Zones zones(2);
	const std::array<std::array<double, 2>, 4> offsets {{{0, 0}, {3, 0},
		{0, 3}, {3, 3}}};
	for (int z {}; z < 4; ++z) add_index_zone(zones, 3, 3, 2);
	auto connect = [&](int a, Zones::Face face, int b)
	{
		zones.set_face(a, face, b);
		zones.set_face(b, static_cast<Zones::Face>(face ^ 1), a);
	};
	connect(0, Zones::xmax, 1);
	connect(2, Zones::xmax, 3);
	connect(0, Zones::ymax, 2);
	connect(1, Zones::ymax, 3);

	Background::Background bkg {};
	std::vector<double> times {0.0, 1.0};
	bkg.move_into_times(times);
	bkg.move_into_zones(zones);
	auto ne {make_field(bkg.get_zones(), 2, [&](int z, double x, double y)
		{return linear(x + offsets[z][0], y + offsets[z][1]);})};
	bkg.move_into_ne(ne);

	// Points within half a cell of the corner at global (3, 3)
	const std::vector<std::array<double, 3>> points {{0, 2.8, 2.9},
		{1, 0.2, 2.7}, {2, 2.6, 0.3}, {3, 0.1, 0.1}, {3, 0.0, 0.0}};
	for (const auto& [zd, x, y] : points)
	{
		const int z {static_cast<int>(zd)};
		EXPECT_NEAR(bkg.interp_ne(z, 0.5, x, y, 1.0),
			linear(x + offsets[z][0], y + offsets[z][1]), 1e-12)
			<< "zone " << z << ", (" << x << ", " << y << ")";
	}
}
