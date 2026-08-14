#include <gtest/gtest.h>
#include <cuda_runtime.h>

#include "utilities.cuh"


// Kernel to call trilinear_interpolate
__global__ void trilinear_test_kernel(double* out)
{
    // Define a simple cube from (0,0,0) to (1,1,1)
    double x0 = 0.0, y0 = 0.0, z0 = 0.0;
    double dx = 1.0, dy = 1.0, dz = 1.0;

    // Vertex values: v = x + y + z at each corner
    double v000 = 0.0;  // (0,0,0)
    double v100 = 1.0;  // (1,0,0)
    double v010 = 1.0;  // (0,1,0)
    double v110 = 2.0;  // (1,1,0)
    double v001 = 1.0;  // (0,0,1)
    double v101 = 2.0;  // (1,0,1)
    double v011 = 2.0;  // (0,1,1)
    double v111 = 3.0;  // (1,1,1)

    // Interpolate at the center of the cube (0.5, 0.5, 0.5)
    double x = 0.5, y = 0.5, z = 0.5;

    double result = Utilities::trilinear_interpolate(
        x0, y0, z0, dx, dy, dz,
        v000, v100, v010, v110,
        v001, v101, v011, v111,
        x, y, z);

    // Write result to device memory
    *out = result;
}


TEST(TrilinearInterpolateCUDA, CenterOfCube)
{
    double* d_out;
    double h_out;

    // Allocate device memory
    cudaMalloc(&d_out, sizeof(double));

    // Launch kernel
    trilinear_test_kernel<<<1,1>>>(d_out);
    cudaDeviceSynchronize();

    // Copy result back
    cudaMemcpy(&h_out, d_out, sizeof(double), cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(d_out);

    // Expected value: at (0.5,0.5,0.5), v = x + y + z = 1.5
    EXPECT_NEAR(h_out, 1.5, 1e-12);
}
