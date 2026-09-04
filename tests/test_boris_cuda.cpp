#include <gtest/gtest.h>

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

	// Slab geometry
	//    | cells | range
	// t  |   128 | 0 - 50 us
	// x  |    32 | 0 - 0.05 m
	// y  |    32 | -0.025, 0.025 m
	// z  |     9 | -0.015, 0.015 m
	//
	// E  = 0 V/m
	// BZ = 1 T   
    opts.set_bkg_source("test");
    opts.set_test_opt("gyrate");

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
    slots.set_z(0, 0.0);     // middle of -0.015–0.015 m

    // Initial velocity: purely along +x
    slots.set_vX(0, 10000.0);
    slots.set_vY(0, 0.0);
    slots.set_vZ(0, 0.0);

    // Charge/mass ratio: electron
    slots.set_q(0, 1.0);      // code multiplies by -e
    slots.set_mass(1.0);

    // Time index and coordinate indices
    slots.set_tidx(0, 0);
    slots.set_xidx(0, 16);
    slots.set_yidx(0, 16);
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
    // Uniform B = (0,0,1) Tesla
    // E = 0
    //
    // Boris rotation angle:
    //     theta = (q/m) * B * dt
    //
    // In your code:
    //     q_m = q * -e / m
    //
    // So q=1 → q_m = -e
    //
    // Therefore:
    //     theta = -e * 1 * dt
    //
    // The velocity rotates in the x–y plane:
    //     v_x' = v_x cos(theta) - v_y sin(theta)
    //     v_y' = v_x sin(theta) + v_y cos(theta)
    //
    // Initial v = (10000, 0, 0)
    //
    // Expected:
    //     v_x' = 10000 * cos(theta)
    //     v_y' = 10000 * sin(theta)
    //

    double q_m = -Constants::charge_e;  // q=1, m=1
    double theta = q_m * dt;            // rad

    double vx_expected = 10000 * std::cos(theta);
    double vy_expected = -10000 * std::sin(theta);

    // -----------------------------
    // 5. Compare with computed
    // -----------------------------
    double vx = slots.vX()[0];
    double vy = slots.vY()[0];
    double vz = slots.vZ()[0];

	std::cout << "vx = " << vx << "  expected = " << vx_expected << "\n";
	std::cout << "vy = " << vy << "  expected = " << vy_expected << "\n";
	std::cout << "vz = " << vz << "  expected = " << 0.0 << "\n";

    EXPECT_NEAR(vx, vx_expected, 1e-6);
    EXPECT_NEAR(vy, vy_expected, 1e-6);
    EXPECT_NEAR(vz, 0.0,        1e-12);

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

