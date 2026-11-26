#include <algorithm>
#include <array>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <unistd.h>

#include "boris.h"
#include "multizone_cases.h"

#include "background.h"
#include "background_device.h"
#include "boris.cuh"
#include "read_bkg.h"
#include "slots.h"
#include "slots_device.h"
#include "options.h"
#include "constants.h"
#include <cmath>

TEST(UpdateVelocityGPU, UniformFieldBorisRotation)
{
    // -----------------------------
    // 1. Options + background
    // -----------------------------
    Options::Options opts {};

	// Slab background (see Slab::read_slab)
	//    | cells | range
	// t  |   128 | 0 - 50 us
	// x  |     8 | 0 - 0.05 m
	// y  |     8 | -0.025, 0.025 m
	// z  |     8 | -0.025, 0.025 m
	//
	// E  = 0 V/m
	// BX = 1 T
    opts.set_bkg_source("slab");
    opts.set_slab_bmag(1.0);

    Background::Background bkg = Background::read_bkg(opts);

	// Copy to device
	Background::BackgroundDevice bkg_d {bkg.to_device()};

    // -----------------------------
    // 2. Create a single particle
    // -----------------------------
    Slots::Slots slots(1);

    // Put particle in the center of the domain
    slots.set_x(0, 0.025);   // middle of 0–0.05 m
    slots.set_y(0, 0.0);     // middle of -0.025–0.025 m
    slots.set_z(0, 0.0);     // middle of -0.025–0.025 m

    // Initial velocity: purely along +Y, across the field
    slots.set_vX(0, 0.0);
    slots.set_vY(0, 10000.0);
    slots.set_vZ(0, 0.0);

    // Charge/mass ratio: electron
    slots.set_q(0, 1.0);      // code multiplies by -e
    slots.set_mass(1.0);

    // Time index and coordinate indices
    slots.set_tidx(0, 0);
    slots.set_xidx(0, 4);
    slots.set_yidx(0, 4);
    slots.set_zidx(0, 4);

    // Actual physical time of particle
    slots.set_t(0, bkg.get_times()[0]);

	// Set state to alive
	slots.set_state(0, 0);

	// Copy to device
	Slots::SlotsDevice slots_d {slots.to_device()};

    // -----------------------------
    // 3. Run update
    // -----------------------------
    double dt = 1e-3;  // 1 ms timestep
	Boris::update_velocity_gpu(slots_d, bkg_d, dt);

	// Copy back to host, copies slots_d into slots. Free memory on GPU.
	slots = slots.to_host(slots_d);
	Slots::free_slots(slots_d);
	Background::free_bkg(bkg_d);

    // -----------------------------
    // 4. Analytic Boris rotation
    // -----------------------------
    //
    // Uniform B = (1,0,0) Tesla
    // E = 0
    //
    // Boris rotation angle:
    //     theta = (q/m) * B * dt
    //
    // In the code:
    //     q_m = q * -e / m
    //
    // So q=1 -> q_m = -e
    //
    // Therefore:
    //     theta = -e * 1 * dt
    //
    // The velocity rotates in the Y-Z plane (dv/dt = q_m v x B, and
    // Y x X = -Z):
    //     v_Y' = v_Y cos(theta)
    //     v_Z' = -v_Y sin(theta)
    //
    // Initial v = (0, 10000, 0)
    //

    double q_m = -Constants::charge_e;  // q=1, m=1
    double theta = q_m * dt;            // rad

    double vy_expected = 10000 * std::cos(theta);
    double vz_expected = -10000 * std::sin(theta);

    // -----------------------------
    // 5. Compare with computed
    // -----------------------------
    double vx = slots.vX()[0];
    double vy = slots.vY()[0];
    double vz = slots.vZ()[0];

	std::cout << "vx = " << vx << "  expected = " << 0.0 << "\n";
	std::cout << "vy = " << vy << "  expected = " << vy_expected << "\n";
	std::cout << "vz = " << vz << "  expected = " << vz_expected << "\n";

    EXPECT_NEAR(vx, 0.0,        1e-12);
    EXPECT_NEAR(vy, vy_expected, 1e-6);
    EXPECT_NEAR(vz, vz_expected, 1e-6);

    // -----------------------------
    // 6. Check reciprocal basis transform
    // -----------------------------
    //
    // In uniform Cartesian background, reciprocal basis should be identity.
    //
    EXPECT_NEAR(slots.vx()[0], vx, 1e-6);
    EXPECT_NEAR(slots.vy()[0], vy, 1e-6);
    EXPECT_NEAR(slots.vz()[0], vz, 1e-6);

}


// ---------------------------------------------------------------------
// The GPU push against the CPU's (Boris::update_velocity_cpu) on the same
// particles, spread over every zone of a multi-zone background: the
// stencils reach into ghost cells at the zone faces, and on the HERMES
// grid E is interpolated along the field line and rotated to the
// particle's toroidal angle. The two agree to rounding (the GPU's fused
// multiply-adds).
// ---------------------------------------------------------------------

namespace
{
	void push_matches_cpu(const Background::Background& bkg_in)
	{
		Background::Background bkg {bkg_in};
		Slots::Slots cpu {MultiZone::spread_particles(bkg, 400, 7)};
		Slots::Slots gpu {cpu};
		const Options::Options opts {};
		const double dt {1e-8};

		Boris::update_velocity_cpu(cpu, bkg, opts, dt, nullptr, nullptr);

		Background::BackgroundDevice bkg_d {bkg.to_device()};
		Slots::SlotsDevice slots_d {gpu.to_device()};
		Boris::update_velocity_gpu(slots_d, bkg_d, dt);
		gpu = gpu.to_host(slots_d);
		Slots::free_slots(slots_d);
		Background::free_bkg(bkg_d);

		// Tolerances relative to the size of each velocity: the Cartesian
		// one, and the computational one (whose components have different
		// units, so each is compared to its own largest value)
		double max_c {};
		std::array<double, 3> max_v {};
		for (int i {}; i < cpu.N(); ++i)
		{
			max_c = std::max({max_c, std::abs(cpu.vX()[i]),
				std::abs(cpu.vY()[i]), std::abs(cpu.vZ()[i])});
			max_v[0] = std::max(max_v[0], std::abs(cpu.vx()[i]));
			max_v[1] = std::max(max_v[1], std::abs(cpu.vy()[i]));
			max_v[2] = std::max(max_v[2], std::abs(cpu.vz()[i]));
		}
		ASSERT_GT(max_c, 0.0);
		int bad {};
		for (int i {}; i < cpu.N(); ++i)
		{
			const bool ok {std::abs(gpu.vX()[i] - cpu.vX()[i]) < 1e-10 * max_c
				&& std::abs(gpu.vY()[i] - cpu.vY()[i]) < 1e-10 * max_c
				&& std::abs(gpu.vZ()[i] - cpu.vZ()[i]) < 1e-10 * max_c
				&& std::abs(gpu.vx()[i] - cpu.vx()[i]) < 1e-10 * max_v[0]
				&& std::abs(gpu.vy()[i] - cpu.vy()[i]) < 1e-10 * max_v[1]
				&& std::abs(gpu.vz()[i] - cpu.vz()[i]) < 1e-10 * max_v[2]};
			if (!ok && bad++ < 5)
			{
				ADD_FAILURE() << "particle " << i << " zone "
					<< cpu.zone()[i] << ": GPU v (" << gpu.vX()[i] << ", "
					<< gpu.vY()[i] << ", " << gpu.vZ()[i] << "), CPU ("
					<< cpu.vX()[i] << ", " << cpu.vY()[i] << ", "
					<< cpu.vZ()[i] << ")";
			}
		}
		EXPECT_EQ(bad, 0);
	}
}

TEST(UpdateVelocityGPU, MatchesCPUOnHermesGrid)
{
	const auto dir {std::filesystem::temp_directory_path()
		/ ("flan_test_boris_cuda_hermes_" + std::to_string(getpid()))};
	const Background::Background bkg {MultiZone::hermes_background(dir)};
	std::filesystem::remove_all(dir);

	// Every zone has a zShift table, so E is field aligned throughout
	for (int zone {}; zone < bkg.get_zones().get_nzones(); ++zone)
		ASSERT_TRUE(bkg.get_zones().has_z_shift(zone));
	push_matches_cpu(bkg);
}

TEST(UpdateVelocityGPU, MatchesCPUOnS3XGrid)
{
	const auto dir {std::filesystem::temp_directory_path()
		/ ("flan_test_boris_cuda_s3x_" + std::to_string(getpid()))};
	const Background::Background bkg {MultiZone::s3x_background(dir)};
	std::filesystem::remove_all(dir);
	ASSERT_GT(bkg.get_zones().get_nzones(), 1);
	push_matches_cpu(bkg);
}
