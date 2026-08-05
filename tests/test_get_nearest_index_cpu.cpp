/*
* AI-written so this is a bit overkill but doesn't hurt.
*/

#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include "impurity_transport.h"


TEST(GetNearestIndexCPU, InteriorPoints)
{
    std::vector<double> vec = {0.0, 1.0, 2.0, 3.0};

    EXPECT_EQ(ImpurityTransport::get_nearest_index_cpu(vec, 0.2), 0);
    EXPECT_EQ(ImpurityTransport::get_nearest_index_cpu(vec, 0.8), 1);
    EXPECT_EQ(ImpurityTransport::get_nearest_index_cpu(vec, 1.6), 2);
    EXPECT_EQ(ImpurityTransport::get_nearest_index_cpu(vec, 2.9), 3);
}

TEST(GetNearestIndexCPU, ExactMatches)
{
    std::vector<double> vec = {0.0, 1.0, 2.0, 3.0};

    EXPECT_EQ(ImpurityTransport::get_nearest_index_cpu(vec, 0.0), 0);
    EXPECT_EQ(ImpurityTransport::get_nearest_index_cpu(vec, 1.0), 1);
    EXPECT_EQ(ImpurityTransport::get_nearest_index_cpu(vec, 2.0), 2);
    EXPECT_EQ(ImpurityTransport::get_nearest_index_cpu(vec, 3.0), 3);
}

TEST(GetNearestIndexCPU, LeftBoundary)
{
    std::vector<double> vec = {0.0, 1.0, 2.0, 3.0};

    EXPECT_EQ(ImpurityTransport::get_nearest_index_cpu(vec, -10.0), 0);
    EXPECT_EQ(ImpurityTransport::get_nearest_index_cpu(vec, -0.1), 0);
}

TEST(GetNearestIndexCPU, RightBoundary)
{
    std::vector<double> vec = {0.0, 1.0, 2.0, 3.0};

    EXPECT_EQ(ImpurityTransport::get_nearest_index_cpu(vec, 10.0), 3);
    EXPECT_EQ(ImpurityTransport::get_nearest_index_cpu(vec, 3.1), 3);
}

TEST(GetNearestIndexCPU, TieBreaks)
{
    std::vector<double> vec = {0.0, 1.0, 2.0, 3.0};

    // Exactly halfway between 1.0 and 2.0 → lower_bound points to 2.0
    // fabs(2.0 - 1.5) == fabs(1.0 - 1.5)
    // lower wins because < is strict
    EXPECT_EQ(ImpurityTransport::get_nearest_index_cpu(vec, 1.5), 1);

    // Slightly closer to 1.0
    EXPECT_EQ(ImpurityTransport::get_nearest_index_cpu(vec, 1.49), 1);

    // Slightly closer to 2.0
    EXPECT_EQ(ImpurityTransport::get_nearest_index_cpu(vec, 1.51), 2);
}

