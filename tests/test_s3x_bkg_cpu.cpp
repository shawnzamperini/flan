/*
* Tests for using a SOLEDGE3X background in a run: loading it through
* read_bkg, starting particles in a chosen zone, and absorbing particles
* that enter the wall. They use the small SOLEDGE3X case written by
* s3x_test_files.h, whose zone 3 has chi = 1 in its last psi row.
*/

#include <filesystem>
#include <gtest/gtest.h>
#include <omp.h>
#include <string>
#include <vector>

#include "background.h"
#include "boundary.h"
#include "options.h"
#include "pcg32.h"
#include "read_bkg.h"
#include "read_s3xe.h"
#include "s3x_test_files.h"
#include "slots.h"


namespace
{
	class S3XBkg : public testing::Test
	{
	protected:
		S3XTest::TestMesh m_mesh {};
		// Named after the test rather than the process, since a death test
		// runs the test again in a new process (see SetUp), which then uses
		// the same files
		std::filesystem::path m_dir {std::filesystem::temp_directory_path()
			/ ("flan_test_s3x_bkg_" + std::string(testing::UnitTest::
			GetInstance()->current_test_info()->name()))};

		void SetUp() override
		{
			// Death tests here read whole backgrounds, which use OpenMP. A
			// forked child of a process whose OpenMP threads already ran
			// can deadlock, so death tests run in a fresh process instead.
			GTEST_FLAG_SET(death_test_style, "threadsafe");

			std::filesystem::create_directories(m_dir);
			const int nphi {m_mesh.nphi};
			S3XTest::write_case(m_dir.string(), m_mesh,
				[nphi](const std::string& path)
				{return S3X::read_mesh(path, nphi).zones;});
		}

		void TearDown() override {std::filesystem::remove_all(m_dir);}

		Options::Options options()
		{
			Options::Options opts {};
			opts.set_bkg_source("soledge3x");
			opts.set_s3x_dir(m_dir.string());
			opts.set_s3x_frame_start(0);
			opts.set_s3x_frame_end(1);
			return opts;
		}
	};

	// One RNG per OpenMP thread, the same as main_loop
	std::vector<pcg32> make_rngs()
	{
		std::vector<pcg32> rngs {};
		for (int i {}; i < omp_get_max_threads(); ++i) rngs.emplace_back(5, i);
		return rngs;
	}

	// One particle in zone (0-based) at (x, y, z), at a time inside the
	// background's
	Slots::Slots one_particle(const Background::Background& bkg, int zone,
		double x, double y, double z)
	{
		Slots::Slots slots(1);
		slots.set_t(0, 0.5 * (bkg.get_times()[0] + bkg.get_times()[1]));
		slots.set_x(0, x);
		slots.set_y(0, y);
		slots.set_z(0, z);
		slots.set_zone(0, zone);
		return slots;
	}
}

// bkg_source = "soledge3x" loads the same background read_s3xe does
TEST_F(S3XBkg, ReadBkg)
{
	const auto opts {options()};
	const auto bkg {Background::read_bkg(opts)};
	const auto direct {S3X::read_s3xe(opts)};

	ASSERT_EQ(bkg.get_nzones(), 3);
	EXPECT_EQ(bkg.get_times(), direct.get_times());
	EXPECT_EQ(bkg.get_ne().get_data(), direct.get_ne().get_data());
	EXPECT_EQ(bkg.get_uX().get_data(), direct.get_uX().get_data());
	EXPECT_EQ(bkg.get_chi().get_data(), direct.get_chi().get_data());
	EXPECT_EQ(bkg.get_zones().get_face(0, Zones::xmax), 2);
}

// s3x_dir has to be set
TEST_F(S3XBkg, UnsetDir)
{
	Options::Options opts {};
	opts.set_bkg_source("soledge3x");
	EXPECT_EXIT(Background::read_bkg(opts),
		testing::ExitedWithCode(EXIT_FAILURE), "s3x_dir must be set");
}

// New particles start in zone imp_zone_start_val (numbered from 1), at a
// position in that zone's coordinates
TEST_F(S3XBkg, StartZone)
{
	auto opts {options()};
	opts.set_imp_zone_start_val(3);
	opts.set_imp_xstart_opt("single_value");
	opts.set_imp_xstart_val(2.5);
	opts.set_imp_ystart_opt("single_value");
	opts.set_imp_ystart_val(10.5);
	opts.set_imp_zstart_opt("single_value");
	opts.set_imp_zstart_val(1.0);
	const auto bkg {Background::read_bkg(opts)};

	Slots::Slots slots(4);
	for (int i {}; i < 4; ++i) slots.set_state(i, 1);
	int rem {4};
	int alive {};
	auto rngs {make_rngs()};
	Slots::fill_slots_cpu(slots, rem, alive, bkg, opts, rngs);

	ASSERT_EQ(alive, 4);
	for (int i {}; i < 4; ++i)
	{
		EXPECT_EQ(slots.zone()[i], 2);
		EXPECT_DOUBLE_EQ(slots.x()[i], 2.5);
		EXPECT_DOUBLE_EQ(slots.y()[i], 10.5);
		EXPECT_EQ(slots.xidx()[i], 2);
		EXPECT_EQ(slots.yidx()[i], 10);
		EXPECT_EQ(slots.zidx()[i], 0);
	}
}

// The start zone has to be one of the background's zones
TEST_F(S3XBkg, StartZoneOutOfRange)
{
	auto opts {options()};
	opts.set_imp_zone_start_val(4);
	EXPECT_EXIT(Background::read_bkg(opts),
		testing::ExitedWithCode(EXIT_FAILURE),
		"imp_zone_start_val = 4, but the background's zones are numbered 1 "
		"to 3");

	opts.set_imp_zone_start_val(0);
	EXPECT_EXIT(Background::read_bkg(opts),
		testing::ExitedWithCode(EXIT_FAILURE), "imp_zone_start_val = 0");

	// A single-zone background only has zone 1
	Options::Options test_opts {};
	test_opts.set_bkg_source("test");
	test_opts.set_imp_zone_start_val(2);
	EXPECT_EXIT(Background::read_bkg(test_opts),
		testing::ExitedWithCode(EXIT_FAILURE), "numbered 1 to 1");
}

// A particle in a cell inside the wall (chi = 1) is absorbed, including
// one that just crossed into it from another zone, unless s3x_chi_absorb
// is off
TEST_F(S3XBkg, WallAbsorbs)
{
	auto opts {options()};
	const auto bkg {Background::read_bkg(opts)};
	auto rngs {make_rngs()};

	auto in_wall {one_particle(bkg, 2, 7.5, 10.5, 1.0)};
	Boundary::Fates fates {};
	Boundary::check_bounds_cpu(in_wall, bkg, opts, rngs, fates);
	EXPECT_EQ(in_wall.state()[0], 1);
	EXPECT_EQ(fates.wall, 1);
	EXPECT_EQ(fates.absorbed, 0);

	auto next_to_wall {one_particle(bkg, 2, 6.5, 10.5, 1.0)};
	Boundary::check_bounds_cpu(next_to_wall, bkg, opts, rngs);
	EXPECT_EQ(next_to_wall.state()[0], 0);

	// From zone 1 across its psi+ face into zone 3's last psi row
	auto crossing {one_particle(bkg, 0, 12.0 + 7.5, 10.5, 1.0)};
	Boundary::check_bounds_cpu(crossing, bkg, opts, rngs);
	EXPECT_EQ(crossing.zone()[0], 2);
	EXPECT_EQ(crossing.state()[0], 1);

	opts.set_s3x_chi_absorb("off");
	auto off {one_particle(bkg, 2, 7.5, 10.5, 1.0)};
	Boundary::check_bounds_cpu(off, bkg, opts, rngs);
	EXPECT_EQ(off.state()[0], 0);
}
