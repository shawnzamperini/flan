#include <cuda_runtime.h>
#include <gtest/gtest.h>

#include "background.h"
#include "background_device.h"
#include "device_constants.cuh"
#include "indices.cuh"
#include "read_bkg.h"
#include "options.h"
#include "interpolate.cuh"


__global__
void test_interp_3d_kernel(Background::BackgroundDevice bkg_d, 
	const double* field, double* result, double test_x, double test_y, 
	double test_z)
{
	// Get indices. d_t, d_x, etc. are from device_constants.cuh.
	int tidx {Indices::get_nearest_index_cuda(d_t, bkg_d.tdim, 0.0)};
	int xidx {Indices::get_nearest_cell_index_cuda(d_grid_x, bkg_d.xdim+1, 
		test_x)};
	int yidx {Indices::get_nearest_cell_index_cuda(d_grid_y, bkg_d.ydim+1, 
		test_y)};
	int zidx {Indices::get_nearest_cell_index_cuda(d_grid_z, bkg_d.zdim+1, 
		test_z)};

    auto s = Interpolate::build_stencil_3d(bkg_d, xidx, yidx, zidx);

    result[0] = Interpolate::interpolate_field_3d(field, s, test_x, test_y, 
		test_z);
}

__global__
void test_interp_4d_kernel(Background::BackgroundDevice bkg_d, 
	const double* field, double* result, double test_t, double test_x, 
	double test_y, double test_z)
{
	// Get indices. d_t, d_x, etc. are from device_constants.cuh.
	int tidx {Indices::get_nearest_index_cuda(d_t, bkg_d.tdim, test_t)};
	int xidx {Indices::get_nearest_cell_index_cuda(d_grid_x, bkg_d.xdim+1, 
		test_x)};
	int yidx {Indices::get_nearest_cell_index_cuda(d_grid_y, bkg_d.ydim+1, 
		test_y)};
	int zidx {Indices::get_nearest_cell_index_cuda(d_grid_z, bkg_d.zdim+1, 
		test_z)};

    auto s = Interpolate::build_stencil_4d(bkg_d, tidx, xidx, yidx, zidx);

    result[0] = Interpolate::interpolate_field_4d(field, s, test_t, test_x, 
		test_y, test_z);
}


TEST(Interpolate3D, LinearField)
{
	// Use one of the test backgrounds, doesn't matter which
    Options::Options opts {};

	// Slab geometry
	//    | cells | range
	// t  |   128 | 0 - 50 us
	// x  |    32 | 0 - 0.05 m
	// y  |    32 | -0.025, 0.025 m
	// z  |     9 | -0.015, 0.015 m
	//
	// E  = 0 V/m
	// BZ = 1 T   
    opts.set_bkg_source("test");
    opts.set_test_opt("gyrate");
    Background::Background bkg = Background::read_bkg(opts);

	// We hijack a 3D field (dxdX) and overwrite it with f(x,y,z) = x + 2y + 3z
	std::cout << bkg.get_dxdX().get_data()[0] << "\n";
	for (int i {}; i < bkg.get_dim2(); i++)
	{
		for (int j {}; j < bkg.get_dim3(); j++)
		{
			for (int k {}; k < bkg.get_dim4(); k++)
			{
				int idx {bkg.get_dxdX().calc_index(i,j,k)};
				double x {bkg.get_x()[i]};
				double y {bkg.get_y()[j]};
				double z {bkg.get_z()[k]};

				bkg.get_dxdX().get_data()[idx] = x + 2 * y + 3 * z;
			}
		}
	}
	std::cout << bkg.get_dxdX().get_data()[0] << "\n";
	

	// Copy to device
	Background::BackgroundDevice bkg_d {bkg.to_device()};

	// Interpolates at (x, y, z) = (0.025, 0.01, 0.01) which equals
	// 0.025 + 2*0.01 + 3*0.01 = 0.224
	double test_x {0.023};
	double test_y {0.09};
	double test_z {0.007};
	double* d_result {};
	cudaMalloc(&d_result, sizeof(double));
    test_interp_3d_kernel<<<1,1>>>(bkg_d, bkg_d.dxdX, d_result, test_x, test_y, 
		test_z);

	// Copy back to host
    double result;
    cudaMemcpy(&result, d_result, sizeof(double), cudaMemcpyDeviceToHost);
	std::cout << "result = " << result << '\n';

	Background::free_bkg(bkg_d);

	double exact {test_x + 2 * test_y + 3 * test_z};
	std::cout << "exact = " << exact << '\n';
    EXPECT_NEAR(result, exact, 1e-12);

}


TEST(Interpolate4D, LinearField)
{
	// Use one of the test backgrounds, doesn't matter which
    Options::Options opts {};

	// Slab geometry
	//    | cells | range
	// t  |   128 | 0 - 50 us
	// x  |    32 | 0 - 0.05 m
	// y  |    32 | -0.025, 0.025 m
	// z  |     9 | -0.015, 0.015 m
	//
	// E  = 0 V/m
	// BZ = 1 T   
    opts.set_bkg_source("test");
    opts.set_test_opt("gyrate");
    Background::Background bkg = Background::read_bkg(opts);

	// We hijack a 4D field (ne) and overwrite it with f(t,x,y,z) = t + x + 2y + 3z
	std::cout << bkg.get_dxdX().get_data()[0] << "\n";
	for (int h {}; h < bkg.get_dim1(); h++)
	{
		for (int i {}; i < bkg.get_dim2(); i++)
		{
			for (int j {}; j < bkg.get_dim3(); j++)
			{
				for (int k {}; k < bkg.get_dim4(); k++)
				{
					int idx {bkg.get_ne().calc_index(h,i,j,k)};
					double t {bkg.get_times()[h]};
					double x {bkg.get_x()[i]};
					double y {bkg.get_y()[j]};
					double z {bkg.get_z()[k]};

					bkg.get_ne().get_data()[idx] = t + x + 2 * y + 3 * z;
				}
			}
		}
	}

	// Copy to device
	Background::BackgroundDevice bkg_d {bkg.to_device()};

	// Interpolates at (t, x, y, z) = (0.000020, 025, 0.01, 0.01) which equals
	// 0.000020 + 0.025 + 2*0.01 + 3*0.01 = 0.22402
	double test_t {0.000020};
	double test_x {0.023};
	double test_y {0.09};
	double test_z {0.007};
	double* d_result {};
	cudaMalloc(&d_result, sizeof(double));
    test_interp_4d_kernel<<<1,1>>>(bkg_d, bkg_d.ne, d_result, test_t, test_x, 
		test_y, test_z);

	// Copy back to host
    double result;
    cudaMemcpy(&result, d_result, sizeof(double), cudaMemcpyDeviceToHost);
	std::cout << "result = " << result << '\n';

	Background::free_bkg(bkg_d);

	double exact {test_t + test_x + 2 * test_y + 3 * test_z};
	std::cout << "exact = " << exact << '\n';
    EXPECT_NEAR(result, exact, 1e-12);

}
