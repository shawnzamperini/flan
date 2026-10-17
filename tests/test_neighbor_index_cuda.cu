#include <gtest/gtest.h>
#include <cuda_runtime.h>

#include "utilities.cuh"


// Kernel to call get_neighbor_index_cuda, since it's device-only
__global__ void neighbor_index_test_kernel(double val,
	const double* centers, int idx, int N, int* out)
{
	*out = Utilities::get_neighbor_index_cuda(val, centers, idx, N);
}


// Copy the cell centers to the device, run the kernel and return the
// neighbor index it found
int get_neighbor_index_on_device(double val, const double* centers,
	int idx, int N)
{
	double* centers_d;
	int* out_d;
	int out {-1};

	cudaMalloc(&centers_d, N * sizeof(double));
	cudaMalloc(&out_d, sizeof(int));
	cudaMemcpy(centers_d, centers, N * sizeof(double),
		cudaMemcpyHostToDevice);

	neighbor_index_test_kernel<<<1,1>>>(val, centers_d, idx, N, out_d);
	EXPECT_EQ(cudaDeviceSynchronize(), cudaSuccess);

	cudaMemcpy(&out, out_d, sizeof(int), cudaMemcpyDeviceToHost);

	cudaFree(centers_d);
	cudaFree(out_d);

	return out;
}


TEST(GetNeighborIndexTest, LeftEdgeMovesRight)
{
    double centers[3] = {0.0, 1.0, 2.0};
    double val = -0.5;   // past leftmost cell center
    int idx = 0;

    int neighbor = get_neighbor_index_on_device(val, centers, idx, 3);

    // At left edge, only valid neighbor is idx + 1
    EXPECT_EQ(neighbor, 1);
}

TEST(GetNeighborIndexTest, RightEdgeMovesLeft)
{
    double centers[3] = {0.0, 1.0, 2.0};
    double val = 2.5;   // past rightmost cell center
    int idx = 2;

    int neighbor = get_neighbor_index_on_device(val, centers, idx, 3);

    // At right edge, only valid neighbor is idx - 1
    EXPECT_EQ(neighbor, 1);
}

TEST(GetNeighborIndexTest, InteriorMovesRight)
{
    double centers[3] = {0.0, 1.0, 2.0};
    double val = 1.2;   // right half of cell 1 center
    int idx = 1;

    int neighbor = get_neighbor_index_on_device(val, centers, idx, 3);

    EXPECT_EQ(neighbor, 2);
}

TEST(GetNeighborIndexTest, InteriorMovesLeft)
{
    double centers[3] = {0.0, 1.0, 2.0};
    double val = 0.7;   // left half of cell 1 center
    int idx = 1;

    int neighbor = get_neighbor_index_on_device(val, centers, idx, 3);

    EXPECT_EQ(neighbor, 0);
}
