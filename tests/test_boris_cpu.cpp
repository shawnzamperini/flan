#include <gtest/gtest.h>

#include "background.h"
#include "boris.h"
#include "read_bkg.h"
#include "slots.h"
#include "options.h"
#include "constants.h"
#include <cmath>

TEST(UpdateVelocityCPU, UniformFieldBorisRotation)
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

    // -----------------------------
    // 3. Run update
    // -----------------------------
    double dt = 1e-3;  // 1 ms timestep
    Boris::update_velocity_cpu(slots, bkg, opts, dt);

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

