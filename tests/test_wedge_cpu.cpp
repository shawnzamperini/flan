/*
* Tests for backgrounds that cover a toroidal wedge rather than the whole
* torus: Cartesian vectors in the ghost cells past the periodic z faces
* are rotated to the ghost's own angle, interpolation across the seam
* is continuous, and a particle that wraps in z has its velocity turned
* with it. They use the small SOLEDGE3X case from s3x_test_files.h, once
* as a quarter of the torus and once as the whole torus.
*/

#include <array>
#include <cmath>
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
	class Wedge : public testing::Test
	{
	protected:
		std::filesystem::path m_dir {std::filesystem::temp_directory_path()
			/ ("flan_test_wedge_" + std::string(testing::UnitTest::
			GetInstance()->current_test_info()->name()))};

		void TearDown() override {std::filesystem::remove_all(m_dir);}

		// Writes the test case spanning phi_span and reads it
		Background::Background read(const double phi_span)
		{
			std::filesystem::remove_all(m_dir);
			std::filesystem::create_directories(m_dir);
			S3XTest::TestMesh mesh {};
			mesh.phi_span = phi_span;
			const int nphi {mesh.nphi};
			S3XTest::write_case(m_dir.string(), mesh,
				[nphi](const std::string& path)
				{return S3X::read_mesh(path, nphi).zones;});
			return Background::read_bkg(options());
		}

		Options::Options options() const
		{
			Options::Options opts {};
			opts.set_bkg_source("soledge3x");
			opts.set_s3x_dir(m_dir.string());
			opts.set_s3x_frame_start(0);
			opts.set_s3x_frame_end(1);
			return opts;
		}
	};

	// (X, Y) rotated about Z by ang
	std::array<double, 2> rotate(const double x, const double y,
		const double ang)
	{
		return {std::cos(ang) * x - std::sin(ang) * y,
			std::sin(ang) * x + std::cos(ang) * y};
	}

	// One RNG per OpenMP thread, the same as main_loop
	std::vector<pcg32> make_rngs()
	{
		std::vector<pcg32> rngs {};
		for (int i {}; i < omp_get_max_threads(); ++i) rngs.emplace_back(5, i);
		return rngs;
	}

	using Field = Vectors::ZonedVector<BkgFPType>;
}

// Every zone of a SOLEDGE3X background is toroidal; a quarter of the
// torus is a wedge and the whole torus isn't
TEST_F(Wedge, IsWedge)
{
	const auto quarter {read(0.5 * S3XTest::pi)};
	const auto whole {read(2.0 * S3XTest::pi)};
	for (int zone {}; zone < 3; ++zone)
	{
		EXPECT_TRUE(quarter.get_zones().is_toroidal(zone));
		EXPECT_TRUE(quarter.get_zones().is_wedge(zone));
		EXPECT_TRUE(whole.get_zones().is_toroidal(zone));
		EXPECT_FALSE(whole.get_zones().is_wedge(zone));
	}
	EXPECT_TRUE(quarter.get_zones().any_wedge());
	EXPECT_FALSE(whole.get_zones().any_wedge());
}

// The field is axisymmetric, so each Cartesian vector, turned back to
// angle 0, is the same in every toroidal cell, including the ghost cells
// past both z faces. Without the rotation, a ghost cell would hold the
// other end's vector, a quarter turn off.
TEST_F(Wedge, GhostVectorsRotated)
{
	const auto bkg {read(0.5 * S3XTest::pi)};
	const auto& zones {bkg.get_zones()};
	const double phi_sign {bkg.get_phi_sign()};

	struct Pair {const Field* x; const Field* y; std::string name;};
	const std::vector<Pair> pairs {
		{&bkg.get_X(), &bkg.get_Y(), "cell center"},
		{&bkg.get_bX(), &bkg.get_bY(), "B"},
		{&bkg.get_gradbX(), &bkg.get_gradbY(), "grad B"},
		{&bkg.get_dzdX(), &bkg.get_dzdY(), "grad z"},
		{&bkg.get_dxdX(), &bkg.get_dxdY(), "grad x"},
		{&bkg.get_eX(), &bkg.get_eY(), "E"},
		{&bkg.get_uX(), &bkg.get_uY(), "u"}};

	for (const auto& p : pairs)
	{
		double scale {};
		for (const double v : p.x->get_data())
			scale = std::max(scale, std::abs(v));
		ASSERT_GT(scale, 0.0) << p.name;

		for (int zone {}; zone < 3; ++zone)
		{
			const int nz {zones.get_nz(zone)};
			const int i {zones.get_nx(zone) / 2};
			const int j {zones.get_ny(zone) / 2};
			auto back = [&](int k)
			{
				const double phi {zones.center(zone, 2, k)};
				return rotate((*p.x)(zone, 0, i, j, k),
					(*p.y)(zone, 0, i, j, k), -phi_sign * phi);
			};
			const auto ref {back(0)};
			for (const int k : {-1, nz / 2, nz})
			{
				const auto v {back(k)};
				EXPECT_NEAR(v[0], ref[0], 1e-9 * scale) << p.name << " k "
					<< k;
				EXPECT_NEAR(v[1], ref[1], 1e-9 * scale) << p.name << " k "
					<< k;
			}
		}
	}
}

// On the whole torus nothing is rotated: the ghost cells are plain copies
TEST_F(Wedge, WholeTorusCopies)
{
	const auto bkg {read(2.0 * S3XTest::pi)};
	const auto& zones {bkg.get_zones()};
	const int nz {zones.get_nz(0)};
	for (const Field* f : {&bkg.get_bX(), &bkg.get_X(), &bkg.get_dzdY()})
	{
		EXPECT_EQ((*f)(0, 0, 3, 5, -1), (*f)(0, 0, 3, 5, nz - 1));
		EXPECT_EQ((*f)(0, 0, 3, 5, nz), (*f)(0, 0, 3, 5, 0));
	}
}

// The two ends of the wedge are the same place a quarter turn apart, so
// B interpolated on the zmax face is B on the zmin face turned by the
// wedge's width
TEST_F(Wedge, InterpolateAcrossSeam)
{
	const double span {0.5 * S3XTest::pi};
	const auto bkg {read(span)};
	const auto& zones {bkg.get_zones()};
	const double t {bkg.get_times()[0]};
	const double x {5.3}, y {12.7};
	const double lo {zones.get_min(0, 2)}, hi {zones.get_max(0, 2)};
	const auto B_lo {bkg.interp_vec3(bkg.get_bX(), bkg.get_bY(),
		bkg.get_bZ(), 0, t, x, y, lo)};
	const auto B_hi {bkg.interp_vec3(bkg.get_bX(), bkg.get_bY(),
		bkg.get_bZ(), 0, t, x, y, hi)};
	const auto turned {rotate(B_lo[0], B_lo[1], bkg.get_phi_sign() * span)};
	const double Bm {std::hypot(B_lo[0], B_lo[1], B_lo[2])};
	EXPECT_NEAR(B_hi[0], turned[0], 1e-9 * Bm);
	EXPECT_NEAR(B_hi[1], turned[1], 1e-9 * Bm);
	EXPECT_NEAR(B_hi[2], B_lo[2], 1e-9 * Bm);
}

// A particle that wraps in z around the wedge has its velocity turned by
// the wedge's width, in the direction it moved. vZ is unchanged.
TEST_F(Wedge, WrapTurnsVelocity)
{
	const double span {0.5 * S3XTest::pi};
	const auto bkg {read(span)};
	const auto opts {options()};
	auto rngs {make_rngs()};
	const double phi_sign {bkg.get_phi_sign()};

	for (const double z : {span + 0.01, -0.01})
	{
		Slots::Slots slots(1);
		slots.set_t(0, bkg.get_times()[0]);
		slots.set_x(0, 5.5);
		slots.set_y(0, 10.5);
		slots.set_z(0, z);
		slots.set_zone(0, 0);
		slots.set_vX(0, 300.0);
		slots.set_vY(0, -200.0);
		slots.set_vZ(0, 50.0);
		Boundary::check_bounds_cpu(slots, bkg, opts, rngs);

		ASSERT_EQ(slots.state()[0], 0);
		const double moved {z > span ? -span : span};
		EXPECT_NEAR(slots.z()[0], z + moved, 1e-12);
		const auto v {rotate(300.0, -200.0, phi_sign * moved)};
		EXPECT_NEAR(slots.vX()[0], v[0], 1e-9);
		EXPECT_NEAR(slots.vY()[0], v[1], 1e-9);
		EXPECT_EQ(slots.vZ()[0], 50.0);
	}
}

// On the whole torus a wrap leaves the velocity exactly as it was
TEST_F(Wedge, WholeTorusWrapKeepsVelocity)
{
	const auto bkg {read(2.0 * S3XTest::pi)};
	const auto opts {options()};
	auto rngs {make_rngs()};
	Slots::Slots slots(1);
	slots.set_t(0, bkg.get_times()[0]);
	slots.set_x(0, 5.5);
	slots.set_y(0, 10.5);
	slots.set_z(0, 2.0 * S3XTest::pi + 0.01);
	slots.set_zone(0, 0);
	slots.set_vX(0, 300.0);
	slots.set_vY(0, -200.0);
	Boundary::check_bounds_cpu(slots, bkg, opts, rngs);
	EXPECT_NEAR(slots.z()[0], 0.01, 1e-12);
	EXPECT_EQ(slots.vX()[0], 300.0);
	EXPECT_EQ(slots.vY()[0], -200.0);
}
