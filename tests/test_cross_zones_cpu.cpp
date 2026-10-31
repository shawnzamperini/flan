/*
* Tests for Boundary::check_bounds_cpu on multi-zone backgrounds: particles
* crossing into neighboring zones, periodic faces, absorbing faces and
* being put back on one of several core faces.
*/

#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <omp.h>
#include <string>
#include <utility>
#include <vector>

#include "background.h"
#include "boundary.h"
#include "options.h"
#include "pcg32.h"
#include "read_bkg.h"
#include "slots.h"
#include "zones.h"


namespace
{
	constexpr double two_pi {6.283185307179586};

	// Grid with n cells between lo and hi, returned as (centers, edges)
	std::pair<std::vector<double>, std::vector<double>> uniform(int n,
		double lo, double hi)
	{
		std::vector<double> centers(n), edges(n + 1);
		const double d {(hi - lo) / n};
		for (int i {}; i <= n; ++i) edges[i] = lo + i * d;
		for (int i {}; i < n; ++i) centers[i] = lo + (i + 0.5) * d;
		return {centers, edges};
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

	// Three zones:
	//
	//        +-------+
	//        |   2   |  ny = 5
	//        +-------+
	//   core |   0   | 1 |  ny = 6
	//        +-------+---+
	//         nx = 4  nx = 3
	//
	// Zone 0's xmax is zone 1's xmin, and zone 0's ymax is zone 2's ymin.
	// Zone 0 and zone 2's xmin are core faces, and every other face in x
	// and y is absorbing.
	Zones::Zones make_three_zones()
	{
		Zones::Zones zones {};
		add_index_zone(zones, 4, 6, 2);
		add_index_zone(zones, 3, 6, 2);
		add_index_zone(zones, 4, 5, 2);
		zones.set_face(0, Zones::xmax, 1);
		zones.set_face(1, Zones::xmin, 0);
		zones.set_face(0, Zones::ymax, 2);
		zones.set_face(2, Zones::ymin, 0);
		zones.set_face(0, Zones::xmin, Zones::face_core);
		zones.set_face(2, Zones::xmin, Zones::face_core);
		return zones;
	}

	// Two zones that wrap around into each other in y, like the core ring
	// of a SOLEDGE3X mesh
	Zones::Zones make_ring()
	{
		Zones::Zones zones {};
		add_index_zone(zones, 3, 4, 2);
		add_index_zone(zones, 3, 7, 2);
		zones.set_face(0, Zones::ymax, 1);
		zones.set_face(1, Zones::ymin, 0);
		zones.set_face(1, Zones::ymax, 0);
		zones.set_face(0, Zones::ymin, 1);
		return zones;
	}

	Background::Background make_background(Zones::Zones zones)
	{
		Background::Background bkg {};
		std::vector<double> times {0.0, 1.0};
		bkg.move_into_times(times);
		bkg.move_into_zones(zones);
		return bkg;
	}

	// One RNG per OpenMP thread, the same as main_loop
	std::vector<pcg32> make_rngs()
	{
		std::vector<pcg32> rngs {};
		for (int i {}; i < omp_get_max_threads(); ++i) rngs.emplace_back(7, i);
		return rngs;
	}

	// A single particle at (x, y, z) in zone, at t = 0.5
	Slots::Slots one_particle(int zone, double x, double y, double z = 1.0)
	{
		Slots::Slots slots(1);
		slots.set_t(0, 0.5);
		slots.set_x(0, x);
		slots.set_y(0, y);
		slots.set_z(0, z);
		slots.set_zone(0, zone);
		return slots;
	}

	// Runs check_bounds_cpu on the slots with the given buffers
	void check_bounds(Slots::Slots& slots, const Background::Background& bkg,
		double buffer = 0.0)
	{
		Options::Options opts {};
		opts.set_imp_xbound_buffer(buffer);
		opts.set_imp_ybound_buffer(buffer);
		opts.set_imp_zbound_buffer(buffer);
		auto rngs {make_rngs()};
		Boundary::check_bounds_cpu(slots, bkg, opts, rngs);
	}
}

// Crossing zone 0's xmax puts the particle at the start of zone 1 in x, and
// leaves y and z alone
TEST(CrossZonesCPU, CrossXmax)
{
	const auto bkg {make_background(make_three_zones())};
	auto slots {one_particle(0, 4.3, 2.5)};
	check_bounds(slots, bkg);

	EXPECT_EQ(slots.state()[0], 0);
	EXPECT_EQ(slots.zone()[0], 1);
	EXPECT_NEAR(slots.x()[0], 0.3, 1e-12);
	EXPECT_DOUBLE_EQ(slots.y()[0], 2.5);
	EXPECT_DOUBLE_EQ(slots.z()[0], 1.0);
}

// And back the other way
TEST(CrossZonesCPU, CrossXmin)
{
	const auto bkg {make_background(make_three_zones())};
	auto slots {one_particle(1, -0.2, 2.5)};
	check_bounds(slots, bkg);

	EXPECT_EQ(slots.zone()[0], 0);
	EXPECT_NEAR(slots.x()[0], 3.8, 1e-12);
}

// Crossing in y
TEST(CrossZonesCPU, CrossYmax)
{
	const auto bkg {make_background(make_three_zones())};
	auto slots {one_particle(0, 1.5, 6.4)};
	check_bounds(slots, bkg);

	EXPECT_EQ(slots.zone()[0], 2);
	EXPECT_DOUBLE_EQ(slots.x()[0], 1.5);
	EXPECT_NEAR(slots.y()[0], 0.4, 1e-12);
}

// A particle past two faces crosses both in one step: out of zone 1
// through xmin into zone 0, then out of zone 0 through ymax into zone 2
TEST(CrossZonesCPU, CornerCrossesTwoFaces)
{
	const auto bkg {make_background(make_three_zones())};
	auto slots {one_particle(1, -0.5, 6.2)};
	check_bounds(slots, bkg);

	EXPECT_EQ(slots.state()[0], 0);
	EXPECT_EQ(slots.zone()[0], 2);
	EXPECT_NEAR(slots.x()[0], 3.5, 1e-12);
	EXPECT_NEAR(slots.y()[0], 0.2, 1e-12);
}

// The z faces connect each zone to itself
TEST(CrossZonesCPU, PeriodicZ)
{
	const auto bkg {make_background(make_three_zones())};
	auto slots {one_particle(1, 1.5, 2.5, two_pi + 0.1)};
	check_bounds(slots, bkg);

	EXPECT_EQ(slots.zone()[0], 1);
	EXPECT_NEAR(slots.z()[0], 0.1, 1e-12);
}

// Absorbing faces kill the particle, counting the buffer. A face shared
// with another zone has no buffer.
TEST(CrossZonesCPU, AbsorbingAndBuffers)
{
	const auto bkg {make_background(make_three_zones())};

	auto past {one_particle(1, 3.2, 2.5)};
	check_bounds(past, bkg);
	EXPECT_EQ(past.state()[0], 1);

	auto in_buffer {one_particle(1, 2.8, 2.5)};
	check_bounds(in_buffer, bkg, 0.3);
	EXPECT_EQ(in_buffer.state()[0], 1);

	auto shared_face {one_particle(0, 3.8, 2.5)};
	check_bounds(shared_face, bkg, 0.5);
	EXPECT_EQ(shared_face.state()[0], 0);
	EXPECT_EQ(shared_face.zone()[0], 0);
	EXPECT_DOUBLE_EQ(shared_face.x()[0], 3.8);
}

// Zones 0 and 2 both have an xmin core face. A particle that hits one is
// put back on either, weighted by their number of cells (6 * 2 = 12 and
// 5 * 2 = 10), at a random y and z within that zone.
TEST(CrossZonesCPU, CoreAcrossZones)
{
	const auto bkg {make_background(make_three_zones())};
	const int n {20000};
	const double buffer {0.1};

	Slots::Slots slots(n);
	for (int i {}; i < n; ++i)
	{
		slots.set_t(i, 0.5);
		slots.set_x(i, -0.05);
		slots.set_y(i, 2.5);
		slots.set_z(i, 1.0);
		slots.set_zone(i, 0);
	}
	check_bounds(slots, bkg, buffer);

	int in_zone2 {};
	for (int i {}; i < n; ++i)
	{
		const int zone {slots.zone()[i]};
		ASSERT_TRUE(zone == 0 || zone == 2) << zone;
		EXPECT_EQ(slots.state()[i], 0);
		EXPECT_DOUBLE_EQ(slots.x()[i], buffer);
		EXPECT_GE(slots.y()[i], buffer);
		EXPECT_LE(slots.y()[i], bkg.get_y_max(zone) - buffer);
		in_zone2 += (zone == 2);
	}

	// Binomial standard deviation is about 70
	EXPECT_NEAR(in_zone2, n * 10.0 / 22.0, 400);
}

// A particle wraps around a ring of zones in either direction
TEST(CrossZonesCPU, Ring)
{
	const auto bkg {make_background(make_ring())};

	auto up {one_particle(1, 1.5, 7.4)};
	check_bounds(up, bkg);
	EXPECT_EQ(up.zone()[0], 0);
	EXPECT_NEAR(up.y()[0], 0.4, 1e-12);

	auto down {one_particle(0, 1.5, -0.4)};
	check_bounds(down, bkg);
	EXPECT_EQ(down.zone()[0], 1);
	EXPECT_NEAR(down.y()[0], 6.6, 1e-12);
}

// A particle far outside the zones keeps crossing until it hits the limit
// and is absorbed instead of looping forever
TEST(CrossZonesCPU, RunawayIsAbsorbed)
{
	const auto bkg {make_background(make_ring())};
	auto slots {one_particle(0, 1.5, 1000.0)};
	check_bounds(slots, bkg);
	EXPECT_EQ(slots.state()[0], 1);
}

// For a single zone, a core face does exactly what core_bc does, drawing
// the same random numbers
TEST(CrossZonesCPU, SingleZoneCoreMatchesCoreBC)
{
	auto [x, gx] = uniform(10, -0.1, 0.1);
	auto [y, gy] = uniform(8, 0.0, 2.0);
	auto [z, gz] = uniform(6, -1.0, 1.0);
	Options::Options opts {};
	opts.set_min_xbound_type("core");
	opts.set_imp_xbound_buffer(0.01);
	opts.set_imp_ybound_buffer(0.02);
	opts.set_imp_zbound_buffer(0.03);
	auto zones {Background::single_zone(x, y, z, gx, gy, gz, opts)};
	const auto bkg {make_background(std::move(zones))};

	Slots::Slots slots(1);
	slots.set_t(0, 0.5);
	slots.set_x(0, -0.095);
	slots.set_y(0, 1.0);
	slots.set_z(0, 0.0);
	std::vector<pcg32> rngs {};
	for (int i {}; i < omp_get_max_threads(); ++i) rngs.emplace_back(3, 0);
	Boundary::check_bounds_cpu(slots, bkg, opts, rngs);

	double a {-0.095}, b {1.0}, c {0.0};
	pcg32 rng(3, 0);
	Boundary::core_bc(a, b, c, -0.1, 0.1, 0.01, 0.0, 2.0, 0.02, -1.0, 1.0,
		0.03, false, rng);

	EXPECT_EQ(slots.state()[0], 0);
	EXPECT_EQ(slots.zone()[0], 0);
	EXPECT_DOUBLE_EQ(slots.x()[0], a);
	EXPECT_DOUBLE_EQ(slots.y()[0], b);
	EXPECT_DOUBLE_EQ(slots.z()[0], c);
}

// Crossing zone 2's ymin puts the particle at the top of zone 0
TEST(CrossZonesCPU, CrossYmin)
{
	const auto bkg {make_background(make_three_zones())};
	auto slots {one_particle(2, 1.5, -0.3)};
	check_bounds(slots, bkg);

	EXPECT_EQ(slots.zone()[0], 0);
	EXPECT_DOUBLE_EQ(slots.x()[0], 1.5);
	EXPECT_NEAR(slots.y()[0], 5.7, 1e-12);
}

// A particle exactly on a shared face hasn't crossed it
TEST(CrossZonesCPU, OnSharedFaceStays)
{
	const auto bkg {make_background(make_three_zones())};

	auto on_max {one_particle(0, 4.0, 2.5)};
	check_bounds(on_max, bkg);
	EXPECT_EQ(on_max.zone()[0], 0);
	EXPECT_DOUBLE_EQ(on_max.x()[0], 4.0);

	auto on_min {one_particle(1, 0.0, 2.5)};
	check_bounds(on_min, bkg);
	EXPECT_EQ(on_min.zone()[0], 1);
	EXPECT_DOUBLE_EQ(on_min.x()[0], 0.0);
}

// After crossing, the new zone's faces are checked too. Here the particle
// lands past (or within the buffer of) zone 1's absorbing xmax.
TEST(CrossZonesCPU, AbsorbedAfterCrossing)
{
	const auto bkg {make_background(make_three_zones())};

	auto past {one_particle(0, 7.5, 2.5)};
	check_bounds(past, bkg);
	EXPECT_EQ(past.state()[0], 1);
	EXPECT_EQ(past.zone()[0], 1);

	auto in_buffer {one_particle(0, 6.9, 2.5)};
	check_bounds(in_buffer, bkg, 0.2);
	EXPECT_EQ(in_buffer.state()[0], 1);
}

// Crossing two zones in the same direction in one step: all the way
// through zone 0 of the ring and back into zone 1
TEST(CrossZonesCPU, TwoCrossingsSameDirection)
{
	const auto bkg {make_background(make_ring())};
	auto slots {one_particle(1, 1.5, 7.0 + 4.0 + 0.5)};
	check_bounds(slots, bkg);

	EXPECT_EQ(slots.state()[0], 0);
	EXPECT_EQ(slots.zone()[0], 1);
	EXPECT_NEAR(slots.y()[0], 0.5, 1e-12);
}

// Dead particles are left alone, even if they're outside their zone
TEST(CrossZonesCPU, DeadParticlesUntouched)
{
	const auto bkg {make_background(make_three_zones())};
	auto slots {one_particle(0, 4.3, -0.5)};
	slots.set_state(0, 1);
	check_bounds(slots, bkg);

	EXPECT_EQ(slots.state()[0], 1);
	EXPECT_EQ(slots.zone()[0], 0);
	EXPECT_DOUBLE_EQ(slots.x()[0], 4.3);
	EXPECT_DOUBLE_EQ(slots.y()[0], -0.5);
}

// The time boundary is checked before the faces. An absorbing one kills
// the particle without moving it, and a periodic one wraps t.
TEST(CrossZonesCPU, TimeBoundary)
{
	const auto bkg {make_background(make_three_zones())};
	auto rngs {make_rngs()};

	Options::Options opts {};
	auto late {one_particle(0, 4.3, 2.5)};
	late.set_t(0, 1.2);
	Boundary::check_bounds_cpu(late, bkg, opts, rngs);
	EXPECT_EQ(late.state()[0], 1);
	EXPECT_EQ(late.zone()[0], 0);
	EXPECT_DOUBLE_EQ(late.x()[0], 4.3);

	opts.set_tbound_type("periodic");
	auto wrapped {one_particle(0, 4.3, 2.5)};
	wrapped.set_t(0, 1.2);
	Boundary::check_bounds_cpu(wrapped, bkg, opts, rngs);
	EXPECT_EQ(wrapped.state()[0], 0);
	EXPECT_NEAR(wrapped.t()[0], 0.2, 1e-12);
	EXPECT_EQ(wrapped.zone()[0], 1);
}

// With several core faces, a particle that hasn't reached one isn't moved
TEST(CrossZonesCPU, CoreNotHit)
{
	const auto bkg {make_background(make_three_zones())};
	auto slots {one_particle(2, 0.3, 2.5)};
	check_bounds(slots, bkg, 0.1);

	EXPECT_EQ(slots.zone()[0], 2);
	EXPECT_DOUBLE_EQ(slots.x()[0], 0.3);
	EXPECT_DOUBLE_EQ(slots.y()[0], 2.5);
	EXPECT_DOUBLE_EQ(slots.z()[0], 1.0);
}

// A particle hitting the core face of the other core zone (zone 2) is put
// back on either core face too, inside that zone's bounds in y and z
TEST(CrossZonesCPU, CoreFromSecondZone)
{
	const auto bkg {make_background(make_three_zones())};
	const int n {2000};
	Slots::Slots slots(n);
	for (int i {}; i < n; ++i)
	{
		slots.set_t(i, 0.5);
		slots.set_x(i, -0.01);
		slots.set_y(i, 4.5);
		slots.set_z(i, 1.0);
		slots.set_zone(i, 2);
	}
	check_bounds(slots, bkg);

	int in_zone0 {};
	for (int i {}; i < n; ++i)
	{
		const int zone {slots.zone()[i]};
		ASSERT_TRUE(zone == 0 || zone == 2) << zone;
		EXPECT_DOUBLE_EQ(slots.x()[i], 0.0);
		EXPECT_LE(slots.y()[i], bkg.get_y_max(zone));
		EXPECT_GE(slots.z()[i], 0.0);
		EXPECT_LE(slots.z()[i], two_pi);
		in_zone0 += (zone == 0);
	}
	EXPECT_GT(in_zone0, 0);
	EXPECT_LT(in_zone0, n);
}

// -------------------------------------------------------------------------
// Single zone: same as before zones
// -------------------------------------------------------------------------

namespace
{
	// The spatial part of check_bounds_cpu as it was before zones (commit
	// d6077ba), for one particle: each bound in the order xmin, xmax,
	// ymin, ymax, zmin, zmax, with the type from the *bound_type options
	// and the outer cell edges as the bounds. Returns the state.
	int old_spatial_bounds(std::array<double, 3>& pos,
		const std::array<int, 6>& types, const std::array<double, 3>& lo,
		const std::array<double, 3>& hi, const std::array<double, 3>& buf,
		pcg32& rng)
	{
		for (int f {}; f < 6; ++f)
		{
			const int a {f / 2};
			const int b {(a + 1) % 3};
			const int c {(a + 2) % 3};
			const bool max_bound {(f & 1) == 1};
			if (types[f] == 0)
			{
				const int state {Boundary::absorbing_bc(pos[a],
					max_bound ? hi[a] : lo[a], max_bound, buf[a])};
				if (state) return state;
			}
			else if (types[f] == 1)
			{
				pos[a] = Boundary::periodic_bc(pos[a], lo[a], hi[a]);
			}
			else if (types[f] == 2)
			{
				Boundary::core_bc(pos[a], pos[b], pos[c], lo[a], hi[a],
					buf[a], lo[b], hi[b], buf[b], lo[c], hi[c], buf[c],
					max_bound, rng);
			}
		}
		return 0;
	}
}

// For every combination of bound types, a single-zone background made
// with single_zone gives exactly what the old bound checks did, including
// the random draws of core boundaries, for particles inside, near and
// outside every bound
TEST(CrossZonesCPU, SingleZoneMatchesOldBounds)
{
	// One thread, so check_bounds_cpu draws from the RNG in particle order
	const int threads {omp_get_max_threads()};
	omp_set_num_threads(1);

	const std::array<std::string, 3> names {"absorbing", "periodic",
		"core"};
	const std::array<double, 3> lo {-0.1, 0.0, -1.0};
	const std::array<double, 3> hi {0.1, 2.0, 1.0};
	const std::array<double, 3> buf {0.01, 0.02, 0.03};
	const int n {60};

	for (int combo {}; combo < 729; ++combo)
	{
		std::array<int, 6> types {};
		int rest {combo};
		for (auto& type : types) {type = rest % 3; rest /= 3;}

		Options::Options opts {};
		opts.set_min_xbound_type(names[types[0]]);
		opts.set_max_xbound_type(names[types[1]]);
		opts.set_min_ybound_type(names[types[2]]);
		opts.set_max_ybound_type(names[types[3]]);
		opts.set_min_zbound_type(names[types[4]]);
		opts.set_max_zbound_type(names[types[5]]);
		opts.set_imp_xbound_buffer(buf[0]);
		opts.set_imp_ybound_buffer(buf[1]);
		opts.set_imp_zbound_buffer(buf[2]);

		auto [x, gx] = uniform(10, lo[0], hi[0]);
		auto [y, gy] = uniform(8, lo[1], hi[1]);
		auto [z, gz] = uniform(6, lo[2], hi[2]);
		auto zones {Background::single_zone(x, y, z, gx, gy, gz, opts)};
		const auto bkg {make_background(std::move(zones))};

		// Positions spread from 30% of the domain below each bound to 30%
		// above it
		pcg32 pos_rng(combo, 1);
		Slots::Slots slots(n);
		std::vector<std::array<double, 3>> want(n);
		for (int i {}; i < n; ++i)
		{
			for (int d {}; d < 3; ++d)
			{
				const double L {hi[d] - lo[d]};
				want[i][d] = lo[d] - 0.3 * L + 1.6 * L * pos_rng.next_double();
			}
			slots.set_t(i, 0.5);
			slots.set_x(i, want[i][0]);
			slots.set_y(i, want[i][1]);
			slots.set_z(i, want[i][2]);
		}

		std::vector<pcg32> rngs {pcg32(combo, 2)};
		Boundary::check_bounds_cpu(slots, bkg, opts, rngs);

		pcg32 rng(combo, 2);
		for (int i {}; i < n; ++i)
		{
			const int state {old_spatial_bounds(want[i], types, lo, hi, buf,
				rng)};
			ASSERT_EQ(slots.state()[i], state) << "combo " << combo
				<< ", particle " << i;
			EXPECT_EQ(slots.zone()[i], 0);
			EXPECT_EQ(slots.x()[i], want[i][0]) << combo << " " << i;
			EXPECT_EQ(slots.y()[i], want[i][1]) << combo << " " << i;
			EXPECT_EQ(slots.z()[i], want[i][2]) << combo << " " << i;
		}
	}

	omp_set_num_threads(threads);
}

// What happens to particles at the boundaries is counted: crossings into
// other zones (two for the corner), absorption at faces, the core and the
// end of the background time, and runaways
TEST(CrossZonesCPU, FatesCounted)
{
	const auto bkg {make_background(make_three_zones())};
	Slots::Slots slots(6);
	const std::vector<std::array<double, 4>> particles {
		{0, 4.3, 2.5, 0.5},		// Crosses into zone 1
		{1, -0.5, 6.2, 0.5},	// Corner, crosses twice into zone 2
		{1, 3.2, 2.5, 0.5},		// Absorbed at zone 1's xmax
		{0, -0.05, 2.5, 0.5},	// Put back in at the core
		{0, 1.5, 2.5, 1.2},		// Past the end of the background time
		{0, 1.5, 2.5, 0.5}};	// Nothing happens
	for (int i {}; i < 6; ++i)
	{
		slots.set_zone(i, static_cast<int>(particles[i][0]));
		slots.set_x(i, particles[i][1]);
		slots.set_y(i, particles[i][2]);
		slots.set_t(i, particles[i][3]);
		slots.set_z(i, 1.0);
	}

	Options::Options opts {};
	auto rngs {make_rngs()};
	Boundary::Fates fates {};
	Boundary::check_bounds_cpu(slots, bkg, opts, rngs, fates);
	EXPECT_EQ(fates.crossings, 3);
	EXPECT_EQ(fates.absorbed, 1);
	EXPECT_EQ(fates.core, 1);
	EXPECT_EQ(fates.time_end, 1);
	EXPECT_EQ(fates.wall, 0);
	EXPECT_EQ(fates.runaway, 0);

	// Counts add up over calls
	const auto ring {make_background(make_ring())};
	auto runaway {one_particle(0, 1.5, 1000.0)};
	Boundary::check_bounds_cpu(runaway, ring, opts, rngs, fates);
	EXPECT_EQ(fates.runaway, 1);
	EXPECT_EQ(fates.crossings, 3 + 8);
}
