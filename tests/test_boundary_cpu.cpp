#include <gtest/gtest.h>
#include "boundary.h"


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
	// their respective bounds (0.0 and 5.0 for both). 
    double a = 12.0, b = 1.0, c = 2.0;
    Boundary::core_bc(a, b, c, 0.0, 10.0, 0.0, 0.0, 5.0, 0.0, 5.0, true);
    EXPECT_LE(a, 10.0);
}
