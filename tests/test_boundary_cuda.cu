#include <gtest/gtest.h>
#include <cuda_runtime.h>

#include "boundary.cuh"
#include "pcg32.h"


// Kernel to call absorbing_bc function
__global__
void absorbing_bc_test_kernel(int* out)
{
    // Case 1: max-bound, value >= bound → dead
    out[0] = Boundary::absorbing_bc_gpu(5.0, 4.0, true);

    // Case 2: max-bound, value < bound → alive
    out[1] = Boundary::absorbing_bc_gpu(3.0, 4.0, true);

    // Case 3: min-bound, value <= bound → dead
    out[2] = Boundary::absorbing_bc_gpu(2.0, 3.0, false);

    // Case 4: min-bound, value > bound → alive
    out[3] = Boundary::absorbing_bc_gpu(5.0, 3.0, false);

    // Case 5: buffer triggers early death
    out[4] = Boundary::absorbing_bc_gpu(3.9, 4.0, true, 0.2);  // 3.9 + 0.2 >= 4.0 → dead
}


// Kernel to call periodic_bc function
__global__
void periodic_bc_test_kernel(double* out)
{
    const double amin = 0.0;
    const double amax = 10.0;

    // Case 1: inside bounds → unchanged
    out[0] = Boundary::periodic_bc_gpu(5.0, amin, amax);

    // Case 2: above max → wrap high
    // 12 > 10 → wrap to 12 - (10 - 0) = 2
    out[1] = Boundary::periodic_bc_gpu(12.0, amin, amax);

    // Case 3: below min → wrap low
    // -3 < 0 → wrap to -3 + (10 - 0) = 7
    out[2] = Boundary::periodic_bc_gpu(-3.0, amin, amax);
}


__global__
void core_bc_test_kernel(double* out_a, double* out_b, double* out_c)
{
    // --- Test parameters ---
    const double a_min = 0.0;
    const double a_max = 10.0;
    const double a_buffer = 1.0;

    const double b_min = -5.0;
    const double b_max =  5.0;
    const double b_buffer = 0.5;

    const double c_min = 100.0;
    const double c_max = 200.0;
    const double c_buffer = 2.0;

    // RNG with deterministic seed
    pcg32 rng(12345u, 67890u);

    // -------------------------------
    // Case 1: HIT max-bound condition
    // -------------------------------
    double a1 = 10.5;   // beyond a_max
    double b1 = 0.0;
    double c1 = 0.0;

    Boundary::core_bc_gpu(a1, b1, c1,
                a_min, a_max, a_buffer,
                b_min, b_max, b_buffer,
                c_min, c_max, c_buffer,
                /*max_bound=*/true,
                rng);

    out_a[0] = a1;
    out_b[0] = b1;
    out_c[0] = c1;

    // -------------------------------
    // Case 2: NO HIT (inside bounds)
    // -------------------------------
    double a2 = 5.0;   // inside
    double b2 = 1.0;
    double c2 = 150.0;

    Boundary::core_bc_gpu(a2, b2, c2,
                a_min, a_max, a_buffer,
                b_min, b_max, b_buffer,
                c_min, c_max, c_buffer,
                /*max_bound=*/true,
                rng);

    out_a[1] = a2;
    out_b[1] = b2;
    out_c[1] = c2;
}


// Test asborbing boundary kills particle
TEST(Boundary, Absorbing)
{
    int* d_out;
    int h_out[5];

    cudaMalloc(&d_out, 5 * sizeof(int));

    absorbing_bc_test_kernel<<<1,1>>>(d_out);
    cudaDeviceSynchronize();

    cudaMemcpy(h_out, d_out, 5 * sizeof(int), cudaMemcpyDeviceToHost);
    cudaFree(d_out);

    EXPECT_EQ(h_out[0], 1);  // dead
    EXPECT_EQ(h_out[1], 0);  // alive
    EXPECT_EQ(h_out[2], 1);  // dead
    EXPECT_EQ(h_out[3], 0);  // alive
    EXPECT_EQ(h_out[4], 1);  // dead due to buffer
}

TEST(PeriodicBCGPU, BasicCases)
{
    double* d_out;
    double h_out[3];

    cudaMalloc(&d_out, 3 * sizeof(double));

    periodic_bc_test_kernel<<<1,1>>>(d_out);
    cudaDeviceSynchronize();

    cudaMemcpy(h_out, d_out, 3 * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_out);

    EXPECT_DOUBLE_EQ(h_out[0], 5.0);  // unchanged
    EXPECT_DOUBLE_EQ(h_out[1], 2.0);  // wrapped high
    EXPECT_DOUBLE_EQ(h_out[2], 7.0);  // wrapped low
}

TEST(CoreBCGPU, BasicCases)
{
    double *d_a, *d_b, *d_c;
    double h_a[2], h_b[2], h_c[2];

    cudaMalloc(&d_a, 2 * sizeof(double));
    cudaMalloc(&d_b, 2 * sizeof(double));
    cudaMalloc(&d_c, 2 * sizeof(double));

    core_bc_test_kernel<<<1,1>>>(d_a, d_b, d_c);
    cudaDeviceSynchronize();

    cudaMemcpy(h_a, d_a, 2 * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_b, d_b, 2 * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_c, d_c, 2 * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);

    // -------------------------------
    // Case 1: HIT max-bound condition
    // -------------------------------
    // a should clamp to a_max - a_buffer = 10 - 1 = 9
    EXPECT_DOUBLE_EQ(h_a[0], 9.0);

    // b should be inside [b_min + b_buffer, b_max - b_buffer] = [-4.5, 4.5]
    EXPECT_GE(h_b[0], -4.5);
    EXPECT_LE(h_b[0],  4.5);

    // c should be inside [c_min + c_buffer, c_max - c_buffer] = [102, 198]
    EXPECT_GE(h_c[0], 102.0);
    EXPECT_LE(h_c[0], 198.0);

    // -------------------------------
    // Case 2: NO HIT
    // -------------------------------
    EXPECT_DOUBLE_EQ(h_a[1], 5.0);
    EXPECT_DOUBLE_EQ(h_b[1], 1.0);
    EXPECT_DOUBLE_EQ(h_c[1], 150.0);
}

