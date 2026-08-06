/*
* AI-written so a bit overkill but doesn't hurt
*/

#include <gtest/gtest.h>
#include <vector>

#include "impurity_transport.h"


// Test interior points
TEST(GetNearestCellIndexCPU, InteriorPoints)
{
    // grid_edges: 0, 1, 2, 3
    // cell centers: 0, 1, 2
    std::vector<double> edges = {0.0, 1.0, 2.0, 3.0};

    EXPECT_EQ(ImpurityTransport::get_nearest_cell_index_cpu(edges, 0.2), 0);  // between 0 and 1
    EXPECT_EQ(ImpurityTransport::get_nearest_cell_index_cpu(edges, 1.4), 1);  // between 1 and 2
    EXPECT_EQ(ImpurityTransport::get_nearest_cell_index_cpu(edges, 2.7), 2);  // between 2 and 3
}

// Test when the value is right at the edge of the grid
TEST(GetNearestCellIndexCPU, ExactEdgeHits)
{
    std::vector<double> edges = {0.0, 1.0, 2.0, 3.0};

    EXPECT_EQ(ImpurityTransport::get_nearest_cell_index_cpu(edges, 1.0), 0);
    EXPECT_EQ(ImpurityTransport::get_nearest_cell_index_cpu(edges, 2.0), 1);
}

// Test when outside the left boundary
TEST(GetNearestCellIndexCPU, LeftBoundary)
{
    std::vector<double> edges = {0.0, 1.0, 2.0, 3.0};

    EXPECT_EQ(ImpurityTransport::get_nearest_cell_index_cpu(edges, -10.0), 0);
    EXPECT_EQ(ImpurityTransport::get_nearest_cell_index_cpu(edges, -0.1), 0);
    EXPECT_EQ(ImpurityTransport::get_nearest_cell_index_cpu(edges, 0.0), 0);
}

// Test when outside the right boundary
TEST(GetNearestCellIndexCPU, RightBoundary)
{
    std::vector<double> edges = {0.0, 1.0, 2.0, 3.0};

    // lower == end() → return index - 2
    // index = 4 → return 2
    EXPECT_EQ(ImpurityTransport::get_nearest_cell_index_cpu(edges, 10.0), 2);
    EXPECT_EQ(ImpurityTransport::get_nearest_cell_index_cpu(edges, 3.1), 2);
    EXPECT_EQ(ImpurityTransport::get_nearest_cell_index_cpu(edges, 3.0), 2);
}

// Just some more tests the AI wrote
TEST(GetNearestCellIndexCPU, VisualSanityCheck)
{
    // Edges: 0   1   2   3
    // Cells: |_0_|_1_|_2_|
    std::vector<double> edges = {0.0, 1.0, 2.0, 3.0};

    // Value between edges[2]=2 and edges[3]=3 → cell index 2
    EXPECT_EQ(ImpurityTransport::get_nearest_cell_index_cpu(edges, 2.5), 2);

    // Value between edges[1]=1 and edges[2]=2 → cell index 1
    EXPECT_EQ(ImpurityTransport::get_nearest_cell_index_cpu(edges, 1.7), 1);

    // Value between edges[0]=0 and edges[1]=1 → cell index 0
    EXPECT_EQ(ImpurityTransport::get_nearest_cell_index_cpu(edges, 0.3), 0);
}

