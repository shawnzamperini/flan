#include <gtest/gtest.h>
#include "utilities.cuh"

TEST(GetNeighborIndexTest, LeftEdgeMovesRight)
{
    double centers[3] = {0.0, 1.0, 2.0};
    double val = -0.5;   // past leftmost cell center
    int idx = 0;

    int neighbor = Utilities::get_neighbor_index_cuda(val, centers, idx, 3);

    // At left edge, only valid neighbor is idx + 1
    EXPECT_EQ(neighbor, 1);
}

TEST(GetNeighborIndexTest, RightEdgeMovesLeft)
{
    double centers[3] = {0.0, 1.0, 2.0};
    double val = 2.5;   // past rightmost cell center
    int idx = 2;

    int neighbor = Utilities::get_neighbor_index_cuda(val, centers, idx, 3);

    // At right edge, only valid neighbor is idx - 1
    EXPECT_EQ(neighbor, 1);
}

TEST(GetNeighborIndexTest, InteriorMovesRight)
{
    double centers[3] = {0.0, 1.0, 2.0};
    double val = 1.2;   // right half of cell 1 center
    int idx = 1;

    int neighbor = Utilities::get_neighbor_index_cuda(val, centers, idx, 3);

    EXPECT_EQ(neighbor, 2);
}

TEST(GetNeighborIndexTest, InteriorMovesLeft)
{
    double centers[3] = {0.0, 1.0, 2.0};
    double val = 0.7;   // left half of cell 1 center
    int idx = 1;

    int neighbor = Utilities::get_neighbor_index_cuda(val, centers, idx, 3);

    EXPECT_EQ(neighbor, 0);
}

