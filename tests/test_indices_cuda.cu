/*
* Tests for get_nearest_index_cuda and get_nearest_cell_index_cuda in
* indices.cuh. The same cases as the CPU tests are checked, and then many
* values on a non-uniform grid are compared against the CPU versions, since
* the two are supposed to give the same results.
*/

#include <gtest/gtest.h>
#include <cuda_runtime.h>
#include <random>
#include <vector>

#include "indices.h"
#include "indices.cuh"


// Kernel to call either get_nearest_index_cuda (cell = false) or
// get_nearest_cell_index_cuda (cell = true) for each value, since they're
// device-only. One thread per value.
__global__ void indices_test_kernel(const double* arr, int N,
	const double* values, int num_values, bool cell, int* out)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= num_values) return;

	if (cell)
		out[i] = Indices::get_nearest_cell_index_cuda(arr, N, values[i]);
	else
		out[i] = Indices::get_nearest_index_cuda(arr, N, values[i]);
}


// Copy arr and values to the device, run the kernel and return the index
// found for each value
std::vector<int> indices_on_device(const std::vector<double>& arr,
	const std::vector<double>& values, bool cell)
{
	const int N = arr.size();
	const int num_values = values.size();

	double* arr_d;
	double* values_d;
	int* out_d;
	cudaMalloc(&arr_d, N * sizeof(double));
	cudaMalloc(&values_d, num_values * sizeof(double));
	cudaMalloc(&out_d, num_values * sizeof(int));
	cudaMemcpy(arr_d, arr.data(), N * sizeof(double),
		cudaMemcpyHostToDevice);
	cudaMemcpy(values_d, values.data(), num_values * sizeof(double),
		cudaMemcpyHostToDevice);

	int blockSize = 256;
	int gridSize = (num_values + blockSize - 1) / blockSize;
	indices_test_kernel<<<gridSize, blockSize>>>(arr_d, N, values_d,
		num_values, cell, out_d);
	EXPECT_EQ(cudaDeviceSynchronize(), cudaSuccess);

	std::vector<int> out(num_values, -1);
	cudaMemcpy(out.data(), out_d, num_values * sizeof(int),
		cudaMemcpyDeviceToHost);

	cudaFree(arr_d);
	cudaFree(values_d);
	cudaFree(out_d);

	return out;
}

std::vector<int> nearest_index(const std::vector<double>& vec,
	const std::vector<double>& values)
{
	return indices_on_device(vec, values, false);
}

std::vector<int> nearest_cell_index(const std::vector<double>& edges,
	const std::vector<double>& values)
{
	return indices_on_device(edges, values, true);
}


// -------------------------------------------------------------------------
// get_nearest_index_cuda
// -------------------------------------------------------------------------

TEST(GetNearestIndexCUDA, InteriorPoints)
{
	std::vector<double> vec = {0.0, 1.0, 2.0, 3.0};
	std::vector<int> expected = {0, 1, 2, 3};

	EXPECT_EQ(nearest_index(vec, {0.2, 0.8, 1.6, 2.9}), expected);
}

TEST(GetNearestIndexCUDA, ExactMatches)
{
	std::vector<double> vec = {0.0, 1.0, 2.0, 3.0};
	std::vector<int> expected = {0, 1, 2, 3};

	EXPECT_EQ(nearest_index(vec, {0.0, 1.0, 2.0, 3.0}), expected);
}

TEST(GetNearestIndexCUDA, OutsideRange)
{
	std::vector<double> vec = {0.0, 1.0, 2.0, 3.0};
	std::vector<int> expected = {0, 0, 3, 3};

	EXPECT_EQ(nearest_index(vec, {-10.0, -0.1, 3.1, 10.0}), expected);
}

TEST(GetNearestIndexCUDA, TieBreaks)
{
	std::vector<double> vec = {0.0, 1.0, 2.0, 3.0};

	// Exactly halfway between 1.0 and 2.0 goes to the lower index
	std::vector<int> expected = {1, 1, 2};

	EXPECT_EQ(nearest_index(vec, {1.5, 1.49, 1.51}), expected);
}

TEST(GetNearestIndexCUDA, SingleValue)
{
	// One background time, e.g. a steady-state background
	std::vector<double> vec = {2.0};
	std::vector<int> expected = {0, 0, 0};

	EXPECT_EQ(nearest_index(vec, {-1.0, 2.0, 5.0}), expected);
}


// -------------------------------------------------------------------------
// get_nearest_cell_index_cuda
// -------------------------------------------------------------------------

TEST(GetNearestCellIndexCUDA, InteriorPoints)
{
	// Edges: 0   1   2   3
	// Cells: |_0_|_1_|_2_|
	std::vector<double> edges = {0.0, 1.0, 2.0, 3.0};
	std::vector<int> expected = {0, 0, 1, 1, 2, 2};

	EXPECT_EQ(nearest_cell_index(edges, {0.2, 0.3, 1.4, 1.7, 2.5, 2.7}),
		expected);
}

TEST(GetNearestCellIndexCUDA, ExactEdgeHits)
{
	// A value on an edge between two cells goes in the lower cell
	std::vector<double> edges = {0.0, 1.0, 2.0, 3.0};
	std::vector<int> expected = {0, 1};

	EXPECT_EQ(nearest_cell_index(edges, {1.0, 2.0}), expected);
}

TEST(GetNearestCellIndexCUDA, LeftBoundary)
{
	std::vector<double> edges = {0.0, 1.0, 2.0, 3.0};
	std::vector<int> expected = {0, 0, 0};

	EXPECT_EQ(nearest_cell_index(edges, {-10.0, -0.1, 0.0}), expected);
}

TEST(GetNearestCellIndexCUDA, RightBoundary)
{
	// N = 4 edges, so the last cell is N - 2 = 2
	std::vector<double> edges = {0.0, 1.0, 2.0, 3.0};
	std::vector<int> expected = {2, 2, 2};

	EXPECT_EQ(nearest_cell_index(edges, {10.0, 3.1, 3.0}), expected);
}

TEST(GetNearestCellIndexCUDA, SingleCell)
{
	std::vector<double> edges = {0.0, 1.0};
	std::vector<int> expected = {0, 0, 0, 0};

	EXPECT_EQ(nearest_cell_index(edges, {-1.0, 0.0, 0.5, 2.0}), expected);
}


// -------------------------------------------------------------------------
// Compare against the CPU versions
// -------------------------------------------------------------------------

// Non-uniform grid edges and the matching cell centers, plus a set of test
// values: random ones over and past the grid, and every edge and center
// exactly (the tie cases).
struct CompareCase
{
	std::vector<double> edges {0.0, 0.1, 0.3, 0.35, 1.0, 2.5, 2.6, 4.0,
		7.0, 7.01, 10.0};
	std::vector<double> centers {};
	std::vector<double> values {};

	CompareCase()
	{
		for (std::size_t i {}; i < edges.size() - 1; ++i)
			centers.push_back(0.5 * (edges[i] + edges[i+1]));

		std::mt19937 gen {12345};
		std::uniform_real_distribution<double> dist {-1.0, 11.0};
		for (int i {}; i < 10000; ++i) values.push_back(dist(gen));

		values.insert(values.end(), edges.begin(), edges.end());
		values.insert(values.end(), centers.begin(), centers.end());
	}
};

TEST(GetNearestIndexCUDA, MatchesCPU)
{
	CompareCase c {};
	std::vector<int> gpu {nearest_index(c.centers, c.values)};

	for (std::size_t i {}; i < c.values.size(); ++i)
	{
		EXPECT_EQ(gpu[i], Indices::get_nearest_index_cpu(c.centers,
			c.values[i])) << "value = " << c.values[i];
	}
}

TEST(GetNearestCellIndexCUDA, MatchesCPU)
{
	CompareCase c {};
	std::vector<int> gpu {nearest_cell_index(c.edges, c.values)};

	for (std::size_t i {}; i < c.values.size(); ++i)
	{
		EXPECT_EQ(gpu[i], Indices::get_nearest_cell_index_cpu(c.edges,
			c.values[i])) << "value = " << c.values[i];
	}
}
