/**
* @file test_trilinear_cuda.cu
*
* @brief Tests of Utilities::trilinear_box on the GPU against the same
* function on the host, which is what the CPU's
* Utilities::trilinear_interpolate calls.
*
* (utilities.h can't be included here: nvcc doesn't get MPI's include
* path, and utilities.h includes mpi.h.)
*/
#include <vector>

#include <gtest/gtest.h>
#include <cuda_runtime.h>

#include "pcg32.h"
#include "utilities_hd.h"


/**
* Evaluates trilinear_box for n cases on the device. Each case is 17
* doubles in the order of trilinear_box's arguments.
*
* @param args The cases' arguments, 17 per case.
* @param out One result per case.
* @param n Number of cases.
*/
__global__ void trilinear_kernel(const double* args, double* out, int n)
{
	const int i {static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x)};
	if (i >= n) return;
	const double* a {args + 17 * i};
	out[i] = Utilities::trilinear_box(a[0], a[1], a[2], a[3], a[4], a[5],
		a[6], a[7], a[8], a[9], a[10], a[11], a[12], a[13], a[14], a[15],
		a[16]);
}


/**
* Runs trilinear_kernel on the cases and returns the results.
*
* @param args The cases' arguments, 17 per case.
* @return One result per case.
*/
std::vector<double> run_on_device(const std::vector<double>& args)
{
	const int n {static_cast<int>(args.size() / 17)};
	double* args_d {};
	double* out_d {};
	cudaMalloc(&args_d, args.size() * sizeof(double));
	cudaMalloc(&out_d, n * sizeof(double));
	cudaMemcpy(args_d, args.data(), args.size() * sizeof(double),
		cudaMemcpyHostToDevice);
	trilinear_kernel<<<(n + 127) / 128, 128>>>(args_d, out_d, n);
	std::vector<double> out(n);
	cudaMemcpy(out.data(), out_d, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(args_d);
	cudaFree(out_d);
	return out;
}


// v = x + y + z on the unit cube, at its center
TEST(TrilinearInterpolateCUDA, CenterOfCube)
{
	const std::vector<double> args {0, 0, 0, 1, 1, 1, 0, 1, 1, 2, 1, 2, 2,
		3, 0.5, 0.5, 0.5};
	EXPECT_NEAR(run_on_device(args)[0], 1.5, 1e-12);
}


// Random boxes, both signs of size, points inside and outside the box
// (clamped), and a box of size 0 in x: the GPU agrees with the host to
// rounding. Not exactly: nvcc fuses a + t * (b - a) into one fused
// multiply-add, which rounds once instead of twice.
TEST(TrilinearInterpolateCUDA, MatchesCpu)
{
	pcg32 rng(12345, 1);
	auto u = [&rng](double lo, double hi)
		{return lo + (hi - lo) * rng.next_double();};

	constexpr int n {4096};
	std::vector<double> args {};
	std::vector<double> expected {};
	for (int i {}; i < n; ++i)
	{
		double c[17] {};
		for (int d {}; d < 3; ++d)
		{
			c[d] = u(-10.0, 10.0);
			c[3 + d] = u(0.1, 2.0) * (u(0.0, 1.0) < 0.5 ? -1.0 : 1.0);

			// Up to half a box outside on either side
			c[14 + d] = c[d] + c[3 + d] * u(-0.5, 1.5);
		}
		if (i % 64 == 0)
		{
			c[3] = 0.0;
			c[14] = c[0];
		}
		for (int v {6}; v < 14; ++v) c[v] = u(-100.0, 100.0);

		args.insert(args.end(), c, c + 17);
		expected.push_back(Utilities::trilinear_box(c[0], c[1], c[2], c[3],
			c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13],
			c[14], c[15], c[16]));
	}

	const std::vector<double> got {run_on_device(args)};
	for (int i {}; i < n; ++i)
	{
		// The values are up to 100, so a few ulps of 100
		EXPECT_FALSE(std::isnan(got[i]));
		EXPECT_NEAR(got[i], expected[i], 1e-12);
	}
}
