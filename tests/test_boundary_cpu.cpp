#include <gtest/gtest.h>
#include "boundary.h"
#include "pcg32.h"


// Test asborbing boundary kills particle
TEST(Boundary, AbsorbingKill)
{
	// Value is 10.0. Bound is 5.0 and it's a maximimum bound (true). Should
	// kill particle because 10 > 5, and thus return 1 (dead).
    int result = Boundary::absorbing_bc(10.0, 5.0, true);
    EXPECT_EQ(result, 1);

	// In this case, 3 < 5 so should return 0 (still alive)
    EXPECT_EQ(Boundary::absorbing_bc(3.0, 5.0, true), 0);
}

TEST(Boundary, PeriodicWrap)
{
	// Value is 11.0. Minimum bound is 0.0 and maximum 10.0. Particle should
	// loop past the max bound and overshoot it by 10.0 - 11.0 = 1.0, so it's
	// new position should be 1.0, which is what's returned.
    EXPECT_EQ(Boundary::periodic_bc(11.0, 0.0, 10.0), 1.0);

	// In this case, goes past the minimum bound by 3.0, so it gets looped
	// around to 10.0 - 3.0 = 7.0.
    EXPECT_EQ(Boundary::periodic_bc(-3.0, 0.0, 10.0), 7.0);
}

TEST(Boundary, CoreHit)
{
	// This one is trickier to explain so check comments around the function
	// definition. We are checking a (12.0) against amin and amax (0.0, 10.0). 
	// We are looking at the maximum bound (true). If a goes past the bound, 
	// then it is placed on that bound and b and c are randomized between
	// their respective bounds (0.0 and 5.0 for both). All the buffers are 
	// zero here.
    double a = 12.0, b = 1.0, c = 2.0;
    pcg32 rng(4, 0);
    Boundary::core_bc(a, b, c, 0.0, 10.0, 0.0, 0.0, 5.0, 0.0, 0.0, 5.0, 0.0,
		true, rng);
    EXPECT_LE(a, 10.0);
}

// LCFS boundary. Domain is x in [0, 1] with lcfs_x = 0.5 (x < 0.5 is core,
// assuming x increases outward), y in [0, 4], z in [-2, 2].
namespace
{
	constexpr double lcfs_x {0.5};
	constexpr double y_min {0.0}, y_max {4.0};
	constexpr double z_min {-2.0}, z_max {2.0};
}

TEST(Boundary, LcfsCoreWrapsMax)
{
	// Core particle past z_max wraps to the other side and stays alive. No
	// twist-shift, so y is unchanged.
	double y {1.0}, z {2.5};
	int state = Boundary::lcfs_bc(y, z, 0.2, lcfs_x, y_min, y_max, z_min,
		z_max, 0.0, 0.0, true);
	EXPECT_EQ(state, 0);
	EXPECT_DOUBLE_EQ(z, -1.5);
	EXPECT_DOUBLE_EQ(y, 1.0);
}

TEST(Boundary, LcfsCoreWrapsMin)
{
	double y {1.0}, z {-2.5};
	int state = Boundary::lcfs_bc(y, z, 0.2, lcfs_x, y_min, y_max, z_min,
		z_max, 0.0, 0.0, false);
	EXPECT_EQ(state, 0);
	EXPECT_DOUBLE_EQ(z, 1.5);
	EXPECT_DOUBLE_EQ(y, 1.0);
}

TEST(Boundary, LcfsCoreOnlyChecksItsSide)
{
	// Past z_max, but we're checking the z_min side, so nothing happens.
	double y {1.0}, z {2.5};
	int state = Boundary::lcfs_bc(y, z, 0.2, lcfs_x, y_min, y_max, z_min,
		z_max, 0.0, 0.0, false);
	EXPECT_EQ(state, 0);
	EXPECT_DOUBLE_EQ(z, 2.5);
}

TEST(Boundary, LcfsSolAbsorbs)
{
	// SOL particle past either z bound is killed and not moved.
	double y {1.0}, z {2.5};
	EXPECT_EQ(Boundary::lcfs_bc(y, z, 0.8, lcfs_x, y_min, y_max, z_min,
		z_max, 0.0, 0.0, true), 1);
	EXPECT_DOUBLE_EQ(z, 2.5);

	z = -2.5;
	EXPECT_EQ(Boundary::lcfs_bc(y, z, 0.8, lcfs_x, y_min, y_max, z_min,
		z_max, 0.0, 0.0, false), 1);

	// Inside the domain, alive
	z = 0.0;
	EXPECT_EQ(Boundary::lcfs_bc(y, z, 0.8, lcfs_x, y_min, y_max, z_min,
		z_max, 0.0, 0.0, true), 0);
}

TEST(Boundary, LcfsAtLcfsIsSol)
{
	// Exactly at lcfs_x counts as SOL (absorbed).
	double y {1.0}, z {2.5};
	EXPECT_EQ(Boundary::lcfs_bc(y, z, lcfs_x, lcfs_x, y_min, y_max, z_min,
		z_max, 0.0, 0.0, true), 1);
}

TEST(Boundary, LcfsBufferOnlyInSol)
{
	// z = 1.95 is inside z_max but within a 0.1 buffer. SOL: absorbed.
	double y {1.0}, z {1.95};
	EXPECT_EQ(Boundary::lcfs_bc(y, z, 0.8, lcfs_x, y_min, y_max, z_min,
		z_max, 0.1, 0.0, true), 1);

	// Core: buffer ignored, not wrapped either.
	EXPECT_EQ(Boundary::lcfs_bc(y, z, 0.2, lcfs_x, y_min, y_max, z_min,
		z_max, 0.1, 0.0, true), 0);
	EXPECT_DOUBLE_EQ(z, 1.95);
}

TEST(Boundary, LcfsTwistShiftHook)
{
	// With a nonzero shift, crossing z_max adds it to y (wrapping in y)
	// and crossing z_min subtracts it.
	double y {3.5}, z {2.5};
	Boundary::lcfs_bc(y, z, 0.2, lcfs_x, y_min, y_max, z_min, z_max, 0.0,
		1.0, true);
	EXPECT_DOUBLE_EQ(y, 0.5);  // 3.5 + 1.0 = 4.5 -> wraps to 0.5

	// Going back across z_min undoes it
	z = -2.5;
	Boundary::lcfs_bc(y, z, 0.2, lcfs_x, y_min, y_max, z_min, z_max, 0.0,
		1.0, false);
	EXPECT_DOUBLE_EQ(y, 3.5);

	// The hook currently returns no shift
	EXPECT_DOUBLE_EQ(Boundary::twist_shift_dy(0.2), 0.0);
}
