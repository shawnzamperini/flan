/**
* @file test_interpolate_cuda.cu
*
* @brief Tests of the GPU interpolation stencils (interpolate.cuh) on a
* single-zone slab and on zones with ghost cells, against exact linear
* fields and the CPU's interpolation.
*/
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include <cuda_runtime.h>
#include <gtest/gtest.h>

#include "background.h"
#include "background_device.h"
#include "indices.cuh"
#include "interpolate.cuh"
#include "options.h"
#include "pcg32.h"
#include "read_bkg.h"
#include "zones.h"
#include "zones_hd.h"


namespace
{
	/// A point to interpolate at
	struct Point
	{
		int zone;
		double t, x, y, z;
	};

	/// The 4D field (ne) and the 3D field (dxdX) at a point
	struct Values
	{
		double ne;
		double dxdX;
	};

	/**
	* Interpolates ne and dxdX at each point, finding the time frame and
	* cell the way the main loop does (find_containing_cell_kernel).
	*/
	__global__ void interp_kernel(const Background::BackgroundDevice bkg_d,
		const Point* p, Values* out, int n)
	{
		const int m {static_cast<int>(blockIdx.x * blockDim.x
			+ threadIdx.x)};
		if (m >= n) return;
		const Point q {p[m]};
		const int tidx {Indices::get_nearest_index_cuda(bkg_d.times,
			bkg_d.tdim, q.t)};
		const int xidx {Zones::cell_index(bkg_d.zones, q.zone, 0, q.x)};
		const int yidx {Zones::cell_index(bkg_d.zones, q.zone, 1, q.y)};
		const int zidx {Zones::cell_index(bkg_d.zones, q.zone, 2, q.z)};

		const auto s3 {Interpolate::build_stencil_3d(bkg_d, q.zone, xidx,
			yidx, zidx, q.x, q.y, q.z)};
		const auto s4 {Interpolate::build_stencil_4d(bkg_d, q.zone, tidx,
			xidx, yidx, zidx, q.t, q.x, q.y, q.z)};
		out[m].dxdX = Interpolate::interpolate_field_3d(bkg_d.dxdX, s3, q.x,
			q.y, q.z);
		out[m].ne = Interpolate::interpolate_field_4d(bkg_d.ne, s4, q.t, q.x,
			q.y, q.z);
	}

	/// Copies the background to the GPU, interpolates at the points
	std::vector<Values> run(Background::Background& bkg,
		const std::vector<Point>& points)
	{
		Background::BackgroundDevice bkg_d {bkg.to_device(0)};
		const int n {static_cast<int>(points.size())};
		Point* p_d {};
		Values* out_d {};
		cudaMalloc(&p_d, n * sizeof(Point));
		cudaMalloc(&out_d, n * sizeof(Values));
		cudaMemcpy(p_d, points.data(), n * sizeof(Point),
			cudaMemcpyHostToDevice);
		interp_kernel<<<(n + 127) / 128, 128>>>(bkg_d, p_d, out_d, n);
		std::vector<Values> out(n);
		cudaMemcpy(out.data(), out_d, n * sizeof(Values),
			cudaMemcpyDeviceToHost);
		EXPECT_EQ(cudaGetLastError(), cudaSuccess);
		cudaFree(p_d);
		cudaFree(out_d);
		Background::free_bkg(bkg_d, 0);
		return out;
	}

	constexpr double two_pi {6.283185307179586};

	/// Grid with n cells between lo and hi, as (centers, edges)
	std::pair<std::vector<double>, std::vector<double>> uniform(int n,
		double lo, double hi)
	{
		std::vector<double> edges(n + 1);
		for (int i {}; i <= n; ++i) edges[i] = lo + i * (hi - lo) / n;
		std::vector<double> centers(n);
		for (int i {}; i < n; ++i)
			centers[i] = 0.5 * (edges[i] + edges[i + 1]);
		return {centers, edges};
	}

	/// Where each zone of make_zones starts in a global x
	constexpr std::array<double, 2> x_offsets {0.0, 4.0};

	/// Field linear in the global x and in y, plus 10 per time frame
	double linear(double gx, double y, double t)
	{
		return 1.0 + 2.0 * gx + 3.0 * y + 10.0 * t;
	}

	/**
	* Two zones side by side in x with ghost width 2, index space in x and
	* y, periodic in z, with the other x and y faces absorbing.
	*/
	Zones::Zones make_zones()
	{
		Zones::Zones zones(2);
		for (const int nx : {4, 3})
		{
			auto [x, gx] = uniform(nx, 0.0, nx);
			auto [y, gy] = uniform(3, 0.0, 3.0);
			auto [z, gz] = uniform(4, 0.0, two_pi);
			const int zone {zones.add_zone(x, y, z, gx, gy, gz)};
			zones.set_face(zone, Zones::zmin, zone);
			zones.set_face(zone, Zones::zmax, zone);
		}
		zones.set_face(0, Zones::xmax, 1);
		zones.set_face(1, Zones::xmin, 0);
		zones.validate();
		return zones;
	}

	/// ZonedVector over the zones with dim1 frames, set to f(zone, t, x,
	/// y) at the cell centers and with the ghost cells filled
	template <typename F>
	Vectors::ZonedVector<BkgFPType> field(const Zones::Zones& zones,
		int dim1, F f)
	{
		Vectors::ZonedVector<BkgFPType> v(dim1, zones.get_nx_all(),
			zones.get_ny_all(), zones.get_nz_all(), zones.get_ghost());
		for (int z {}; z < zones.get_nzones(); ++z)
		for (int t {}; t < dim1; ++t)
		for (int i {}; i < zones.get_nx(z); ++i)
		for (int j {}; j < zones.get_ny(z); ++j)
		for (int k {}; k < zones.get_nz(z); ++k)
			v(z, t, i, j, k) = f(z, t, zones.get_x(z)[i], zones.get_y(z)[j]);
		zones.fill_ghosts(v);
		return v;
	}

	/**
	* Background over make_zones with times {0, 1}: ne and dxdX are linear
	* in the global x (ne also in time), and every other field the GPU
	* copies is zero.
	*/
	Background::Background make_background()
	{
		Background::Background bkg {};
		std::vector<double> times {0.0, 1.0};
		bkg.move_into_times(times);
		auto zones {make_zones()};
		bkg.move_into_zones(zones);
		const auto& zs {bkg.get_zones()};

		auto zero4 = [&]() {return field(zs, 2,
			[](int, int, double, double) {return 0.0;});};
		auto zero3 = [&]() {return field(zs, 1,
			[](int, int, double, double) {return 0.0;});};

		auto ne {field(zs, 2, [](int z, int t, double x, double y)
			{return linear(x + x_offsets[z], y, t);})};
		auto dxdX {field(zs, 1, [](int z, int, double x, double y)
			{return linear(x + x_offsets[z], y, 0.0);})};
		bkg.move_into_ne(ne);
		bkg.move_into_dxdX(dxdX);

		{auto f {zero4()}; bkg.move_into_te(f);}
		{auto f {zero4()}; bkg.move_into_ti(f);}
		{auto f {zero4()}; bkg.move_into_bX(f);}
		{auto f {zero4()}; bkg.move_into_bY(f);}
		{auto f {zero4()}; bkg.move_into_bZ(f);}
		{auto f {zero4()}; bkg.move_into_bmag(f);}
		{auto f {zero4()}; bkg.move_into_eX(f);}
		{auto f {zero4()}; bkg.move_into_eY(f);}
		{auto f {zero4()}; bkg.move_into_eZ(f);}
		{auto f {zero4()}; bkg.move_into_emag(f);}
		{auto f {zero4()}; bkg.move_into_uX(f);}
		{auto f {zero4()}; bkg.move_into_uY(f);}
		{auto f {zero4()}; bkg.move_into_uZ(f);}
		{auto f {zero3()}; bkg.move_into_dxdY(f);}
		{auto f {zero3()}; bkg.move_into_dxdZ(f);}
		{auto f {zero3()}; bkg.move_into_dydX(f);}
		{auto f {zero3()}; bkg.move_into_dydY(f);}
		{auto f {zero3()}; bkg.move_into_dydZ(f);}
		{auto f {zero3()}; bkg.move_into_dzdX(f);}
		{auto f {zero3()}; bkg.move_into_dzdY(f);}
		{auto f {zero3()}; bkg.move_into_dzdZ(f);}
		return bkg;
	}
}


// On the slab (one zone, no ghost cells), fields linear in x, y, z and t
// are interpolated exactly
TEST(InterpolateCUDA, SlabLinearField)
{
	// The slab background with its defaults: 128 frames over 0 - 50 us,
	// 8 x 8 x 8 cells over x in 0 - 0.05 m, y and z in -0.025 - 0.025 m
	Options::Options opts {};
	opts.set_bkg_source("slab");
	Background::Background bkg {Background::read_bkg(opts)};

	// dxdX = x + 2y + 3z and ne = t + x + 2y + 3z
	for (int h {}; h < bkg.get_dim1(); ++h)
	for (int i {}; i < bkg.get_dim2(); ++i)
	for (int j {}; j < bkg.get_dim3(); ++j)
	for (int k {}; k < bkg.get_dim4(); ++k)
	{
		const double t {bkg.get_times()[h]};
		const double x {bkg.get_x(0)[i]};
		const double y {bkg.get_y(0)[j]};
		const double z {bkg.get_z(0)[k]};
		bkg.get_ne().get_data()[bkg.get_ne().calc_index(0, h, i, j, k)] =
			t + x + 2 * y + 3 * z;
		if (h == 0)
		{
			bkg.get_dxdX().get_data()[bkg.get_dxdX().calc_index(0, 0, i, j,
				k)] = x + 2 * y + 3 * z;
		}
	}

	const Point p {0, 20e-6, 0.023, 0.01, 0.007};
	const auto got {run(bkg, {p})};
	EXPECT_NEAR(got[0].dxdX, p.x + 2 * p.y + 3 * p.z, 1e-12);
	EXPECT_NEAR(got[0].ne, p.t + p.x + 2 * p.y + 3 * p.z, 1e-12);
}


// With two zones and ghost cells, points near the shared face in either
// zone interpolate through the ghosts: the GPU gets the exact linear field,
// the same as the CPU's interp_ne
TEST(InterpolateCUDA, AcrossZones)
{
	auto bkg {make_background()};
	const auto& zones {bkg.get_zones()};

	pcg32 rng(5, 2);
	std::vector<Point> points {};
	for (int n {}; n < 4000; ++n)
	{
		Point p {};
		p.zone = static_cast<int>(rng.next_double() * 2);
		p.t = rng.next_double();
		p.x = zones.get_min(p.zone, 0) + rng.next_double()
			* (zones.get_max(p.zone, 0) - zones.get_min(p.zone, 0));
		p.y = 0.5 + 2.0 * rng.next_double();
		p.z = two_pi * rng.next_double();
		points.push_back(p);
	}

	const auto got {run(bkg, points)};
	for (std::size_t i {}; i < points.size(); ++i)
	{
		const Point& p {points[i]};
		const double gx {p.x + x_offsets[p.zone]};

		// The outer x faces are absorbing, so the field is held past the
		// first and last centers (global x 0.5 and 6.5)
		const double hx {std::min(std::max(gx, 0.5), 6.5)};
		EXPECT_NEAR(got[i].ne, linear(hx, p.y, p.t), 1e-11) << i;
		EXPECT_NEAR(got[i].dxdX, linear(hx, p.y, 0.0), 1e-11) << i;
		EXPECT_NEAR(got[i].ne, bkg.interp_ne(p.zone, p.t, p.x, p.y, p.z),
			1e-11) << i;
	}
}
