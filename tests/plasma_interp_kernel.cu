/*
* Kernel and launcher for test_plasma_interp_cuda (see
* plasma_interp_kernel.h).
*/

#include <cuda_runtime.h>

#include "plasma_interp_kernel.h"

#include "interpolate.cuh"

namespace PlasmaInterpTest
{
	__global__ void interp_kernel(const Background::BackgroundDevice bkg_d,
		const Point* points, const int n, double* out)
	{
		const int i {static_cast<int>(blockIdx.x * blockDim.x
			+ threadIdx.x)};
		if (i >= n) return;
		const Point p {points[i]};
		const Interpolate::InterpolationStencil s {
			Interpolate::build_stencil_4d(bkg_d, p.zone, p.tidx, p.xidx,
			p.yidx, p.zidx, p.t, p.x, p.y, p.z)};
		double* o {out + static_cast<std::size_t>(i) * num_values};
		o[0] = Interpolate::interp_plasma(bkg_d, bkg_d.ne, s, p.zone, p.t,
			p.x, p.y, p.z);
		o[1] = Interpolate::interp_plasma(bkg_d, bkg_d.te, s, p.zone, p.t,
			p.x, p.y, p.z);
		o[2] = Interpolate::interp_plasma(bkg_d, bkg_d.ti, s, p.zone, p.t,
			p.x, p.y, p.z);
		o[3] = Interpolate::interp_plasma(bkg_d, bkg_d.emag, s, p.zone, p.t,
			p.x, p.y, p.z);
		Interpolate::interp_plasma_vec3(bkg_d, bkg_d.uX, bkg_d.uY, bkg_d.uZ,
			s, p.zone, p.t, p.x, p.y, p.z, o + 4);
	}

	std::vector<double> interp_on_gpu(
		const Background::BackgroundDevice& bkg_d,
		const std::vector<Point>& points)
	{
		const int n {static_cast<int>(points.size())};
		Point* points_d {};
		double* out_d {};
		cudaSetDevice(bkg_d.device_id);
		cudaMalloc(&points_d, n * sizeof(Point));
		cudaMalloc(&out_d, n * num_values * sizeof(double));
		cudaMemcpy(points_d, points.data(), n * sizeof(Point),
			cudaMemcpyHostToDevice);
		interp_kernel<<<(n + 127) / 128, 128>>>(bkg_d, points_d, n, out_d);
		std::vector<double> out (static_cast<std::size_t>(n) * num_values);
		cudaMemcpy(out.data(), out_d, out.size() * sizeof(double),
			cudaMemcpyDeviceToHost);
		cudaFree(points_d);
		cudaFree(out_d);
		return out;
	}
}
