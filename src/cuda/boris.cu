#include <cuda_runtime.h>
#include <cstdio>

#include "background_device.h"
#include "slots_device.h"

#include "device_constants.cuh"
#include "utilities.cuh"


namespace Boris
{

	__global__ void update_velocity_kernel(Slots::SlotsDevice slots_d, 
		const Background::BackgroundDevice bkg_d)
	{
		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;

		// Don't try and access beyond the number of slots (segfault)
		if (i >= slots_d.N) return;

		// Local variables
		int tidx {slots_d.tidx[i]};
		int xidx {slots_d.xidx[i]};
		int yidx {slots_d.yidx[i]};
		int zidx {slots_d.zidx[i]};

		double t {slots_d.t[i]};
		double x {slots_d.x[i]};
		double y {slots_d.y[i]};
		double z {slots_d.z[i]};

		// Get nearest neighbor indices for each direction. These tell us
		// which direction we should interpolate towards, i.e., which
		// rectangle made by the neighboring cell centers our particle
		// is bounded by.
		// d_x is a global constant defined in cuda/device_constants.h.
		const int xidx_neighbor {Utilities::get_neighbor_index_cuda(x, 
			Background::d_x, xidx, bkg_d.xdim)};
		const int yidx_neighbor {Utilities::get_neighbor_index_cuda(y, 
			Background::d_y, yidx, bkg_d.ydim)};
		const int zidx_neighbor {Utilities::get_neighbor_index_cuda(z, 
			Background::d_z, zidx, bkg_d.zdim)};

		// Similarly for t, except we can't use get_neighbor_index since it
		// uses cell center coordinates, and t is defined at each frame (i.e.,
		// not between each frame, that would be nonintuitive). So we use a
		// little SIMD-friendly logic here to assign tidx_neighbor to tidx+1,
		// and tidx-1 if we're in the last time frame.
		// at_end = 1 if tidx == ntimes-1, else 0 
		// If at_end = 0 → i + 1 
		// If at_end = 1 → i - 1 
		const int at_end {(tidx == static_cast<int>(bkg_d.tdim) - 1)}; 
		const int tidx_neighbor {tidx + 1 - 2 * at_end};

		// Time between neighboring t frames
		const double frame_dt {d_t[tidx_neighbor] - d_t[tidx]};

		// Trilinear interpolation to ensure we use a continous B/E in the 
		// algorithm. Without interpolating, the particle will:
		// a) Fail to follow the parallel direction of the field line due to
		//    it being discrete in the code (this is called "numerical
		//    diffusion" and introduces artificial cross-field transport).
		// b) Experience large "kicks" at cell boundaries where B/E is 
		//    discontinuous, which can introduce artificial drifts.
		// 
		// The issue in a) generally will still occur since we are still using
		// an approximation to the field line, but this greatly reduces the
		// numerical diffusion since B is now continous.

		// x, y, z coordinates of two bounding vertices to interpolate 
		// between. Note these are not grid vertices, but rather are formed
		// by cell center coordinates since that's where B/E are assumed
		// to be defined. It is essentially a cell shifted by dx/2, dy/2
		// and dz/2 if that helps.
		const double x0 {d_x[xidx]};
		const double x1 {d_x[xidx_neighbor]};
		const double y0 {d_y[yidx]};
		const double y1 {d_y[yidx_neighbor]};
		const double z0 {d_z[zidx]};
		const double z1 {d_z[zidx_neighbor]};

		// Find a nice way to group and interpolate E and B
		// ...

	}

	void update_velocity_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d, 
		const double dt)
	{

	/*
		#pragma omp parallel for
		for (int i=0; i < slots.N(); ++i)
		{
			// Local variables
			int tidx {slots.tidx()[i]};
			int xidx {slots.xidx()[i]};
			int yidx {slots.yidx()[i]};
			int zidx {slots.zidx()[i]};

			// Get nearest neighbor indices for each direction. These tell us
			// which direction we should interpolate towards, i.e., which
			// rectangle made by the neighboring cell centers our particle
			// is bounded by.
			const int xidx_neighbor {Utilities::get_neighbor_index(slots.x()[i], 
				bkg.get_x(), xidx)};
			const int yidx_neighbor {Utilities::get_neighbor_index(slots.y()[i], 
				bkg.get_y(), yidx)};
			const int zidx_neighbor {Utilities::get_neighbor_index(slots.z()[i], 
				bkg.get_z(), zidx)};

			// Similarly for t, except we can't use get_neighbor_index since it
			// uses cell center coordinates, and t is defined at each frame (i.e.,
			// not between each frame, that would be nonintuitive). So we use a
			// little SIMD-friendly logic here to assign tidx_neighbor to tidx+1,
			// and tidx-1 if we're in the last time frame.
			// at_end = 1 if tidx == ntimes-1, else 0 
			// If at_end = 0 → i + 1 
			// If at_end = 1 → i - 1 
			const int at_end {(tidx == static_cast<int>
				(bkg.get_times().size())-1)}; 
			const int tidx_neighbor {tidx + 1 - 2 * at_end};

			//std::cout << "t: " << tidx << "\t" << tidx_neighbor << "\n";
			//std::cout << "x: " << xidx << "\t" << xidx_neighbor << "\n";
			//std::cout << "y: " << yidx << "\t" << yidx_neighbor << "\n";
			//std::cout << "z: " << zidx << "\t" << zidx_neighbor << "\n";

			// Time between neighboring t frames
			const double frame_dt {bkg.get_times()[tidx_neighbor] 
				- bkg.get_times()[tidx]};

			// Trilinear interpolation to ensure we use a continous B/E in the 
			// algorithm. Without interpolating, the particle will:
			// a) Fail to follow the parallel direction of the field line due to
			//    it being discrete in the code (this is called "numerical
			//    diffusion" and introduces artificial cross-field transport).
			// b) Experience large "kicks" at cell boundaries where B/E is 
			//    discontinuous, which can introduce artificial drifts.
			// 
			// The issue in a) generally will still occur since we are still using
			// an approximation to the field line, but this greatly reduces the
			// numerical diffusion since B is now continous.
			//
			// In the following, we wrap the components into a couple nested arrays
			// so that we can use a single nested loop to iterate through each 
			// B and E component and not have to write the trilinear interpolation
			// loop twice for each vector (or even 6 times, one for each vector
			// component!). It adds a bit of complexity but avoids repeating code.

			// Wrap in convienent array. Need a const reference in reference_wrapper
			// because bkg.get_bX() returns a const ref. Need to match the getter. 
			const std::array<std::reference_wrapper<
				const Vectors::Vector4D<BkgFPType>>, 4> B 
				{bkg.get_bX(), bkg.get_bY(), bkg.get_bZ(), bkg.get_bmag()};
			const std::array<std::reference_wrapper<
				const Vectors::Vector4D<BkgFPType>>, 4> E 
				{bkg.get_eX(), bkg.get_eY(), bkg.get_eZ(), bkg.get_emag()};

			// And put again into single array so we can make an outer loop out
			// of these.
			const std::array outer {&B, &E};

			// Similarly, arrays to hold the output local values
			std::array<double, 4> B_local {};
			std::array<double, 4> E_local {};
			std::array<double, 4>* locals [] = {&B_local, &E_local};

			// x, y, z coordinates of two bounding vertices to interpolate 
			// between. Note these are not grid vertices, but rather are formed
			// by cell center coordinates since that's where B/E are assumed
			// to be defined. It is essentially a cell shifted by dx/2, dy/2
			// and dz/2 if that helps.
			const double x0 {bkg.get_x()[xidx]};
			const double x1 {bkg.get_x()[xidx_neighbor]};
			const double y0 {bkg.get_y()[yidx]};
			const double y1 {bkg.get_y()[yidx_neighbor]};
			const double z0 {bkg.get_z()[zidx]};
			const double z1 {bkg.get_z()[zidx_neighbor]};

			// Loop through each B/E component, performing a trilinear interpolation
			// to get the local value of each. First loop is to choose B or E.
			for (int j = 0; j < 2; ++j) 
			{ 
				// Index out (B or E), dereference and assign to comps. This is 
				// an object the same as B and E above.
				const auto& comps = *outer[j]; 

				// Get reference to output local values. This is a reference to
				// an array of 3 doubles, B_local or E_local.
				auto& out = *locals[j]; 

				// Loop through each component to get interpolated value for each
				for (int k = 0; k < 4; ++k) 
				{
					// Values at each vertex, 8 total because it's a rectangle.
					double v000 {comps[k](tidx, xidx, yidx, zidx)};
					double v100 {comps[k](tidx, xidx_neighbor, yidx, zidx)};
					double v010 {comps[k](tidx, xidx, yidx_neighbor, zidx)};
					double v110 {comps[k](tidx, xidx_neighbor, yidx_neighbor, 
						zidx)};
					double v001 {comps[k](tidx, xidx, yidx, zidx_neighbor)};
					double v101 {comps[k](tidx, xidx_neighbor, yidx, 
						zidx_neighbor)};
					double v011 {comps[k](tidx, xidx, yidx_neighbor, 
						zidx_neighbor)};
					double v111 {comps[k](tidx, xidx_neighbor, yidx_neighbor,	
						zidx_neighbor)};

					// Perform interpolation, storing value in our local array
					const double comp0 = Utilities::trilinear_interpolate(x0, y0, 
						z0, x1, y1, z1, v000, v100, v010, v110, v001, v101, v011, 
						v111, slots.x()[i], slots.y()[i], slots.z()[i]);

					// Now repeat for the previous time index so we can do a simple
					// linear interpolation in time. This is critical to get a
					// decent estimate of the polarization drift, for instance.
					// If we don't do this, the field will be discontinuous each
					// time the particle enters a new time frame, causing 
					// unphysical kicks.
					v000 = comps[k](tidx_neighbor, xidx, yidx, zidx);
					v100 = comps[k](tidx_neighbor, xidx_neighbor, yidx, zidx);
					v010 = comps[k](tidx_neighbor, xidx, yidx_neighbor, zidx);
					v110 = comps[k](tidx_neighbor, xidx_neighbor, 
						yidx_neighbor, zidx);
					v001 = comps[k](tidx_neighbor, xidx, yidx, zidx_neighbor);
					v101 = comps[k](tidx_neighbor, xidx_neighbor, yidx, 
						zidx_neighbor);
					v011 = comps[k](tidx_neighbor, xidx, yidx_neighbor, 
						zidx_neighbor);
					v111 = comps[k](tidx_neighbor, xidx_neighbor, 
						yidx_neighbor, zidx_neighbor);

					// Interpolate again at new time frame
					const double comp1 = Utilities::trilinear_interpolate(x0, y0, 
						z0, x1, y1, z1, v000, v100, v010, v110, v001, v101, v011, 
						v111, slots.x()[i], slots.y()[i], slots.z()[i]);

					// Now linearly interpolate in time (this is just point slope)
					const double slope {(comp1 - comp0) / frame_dt};
					out[k] = slope * (slots.t()[i] - bkg.get_times()[tidx]) 
						+ comp0;
				}  // k loop

				// We interpolated the magnitude as well. This is because there is
				// no guarantee that the interpolated vector will have the same
				// magnitude (in fact it almost certainly shrinks). We are faced
				// with a choice:
				//   a) Leave as-is. B/E are continuous at cell boundaries, but
				//      the magntiude will be off. The size of the gyro-orbit will
				//      be affected as a result.
				//   b) Scale via the interpolated magnitude. B/E will be 
				//      discontinuous at cell boundaries, but definitely not as
				//      discontinous as without interpolation. Gyro-orbit will
				//      be more physically correct.
				// We choose option b), because getting the gyro-orbit correct is
				// central to Flan. Slight discontinuities are not ideal, but the
				// particle spends order of magnitude more time not crossing 
				// a boundary, therefore option b) also minimizes numerical
				// diffusion in that sense as well by ensuring the gyro-orbit is
				// as correct as possible.
				const double out_mag {std::sqrt(out[0] * out[0] + out[1] * out[1] 
					+ out[2] * out[2])};
				double mask {static_cast<double>(out_mag > 0)};
				for (int k {}; k < 3; ++k)
				{
					// out[3] contains the interpolated magnitude. If out_mag
					// is zero, mask will just multiply by zero to force the vector
					// to zero. It's a divide-by-zero guard.
					out[k] = out[k] * out[3] / std::max(out_mag, 1e-30) * mask;
				}  // k loop

			}  // j loop

			// Normal Boris algorithm continues

			// Store in local variables for conveinence
			double q_m {slots.q()[i] * -Constants::charge_e / slots.mass()};
			std::array<double, 3> v {slots.vX()[i], slots.vY()[i], slots.vZ()[i]};

			// t vector
			std::array<double, 3> t {};
			for (int j {}; j < 3; ++j) t[j] = q_m * B_local[j] * 0.5 * dt;

			// Magnitude of t, squared
			double tmag2 {t[0]*t[0] + t[1]*t[1] + t[2]*t[2]};

			// s vector
			std::array<double, 3> s {};
			for (int j {}; j < 3; ++j) s[j] = 2 * t[j] / (1.0 + tmag2);

			// v minus
			std::array<double, 3> vminus {};
			for (int j {}; j < 3; ++j) vminus[j] = v[j] + q_m * E_local[j] * 0.5*dt;

			// v prime
			std::array<double, 3> vprime {};
			std::array<double, 3> vminus_cross_t {Utilities::cross_product(vminus,	
				t)};
			for (int j {}; j < 3; ++j) vprime[j] = vminus[j] + vminus_cross_t[j];

			// v plus
			std::array<double, 3> vplus {};
			std::array<double ,3> vprime_cross_s {Utilities::cross_product(vprime, 
				s)};
			for (int j {}; j < 3; ++j) vplus[j] = vminus[j] + vprime_cross_s[j];

			// v n+1/2
			// Note we are storing particle velocity at the half time steps before
			// the position. I.e., xi = x(ti), vi = v(ti - dt/2). 
			slots.set_vX(i, vplus[0] + q_m * E_local[0] * 0.5 * dt);
			slots.set_vY(i, vplus[1] + q_m * E_local[1] * 0.5 * dt);
			slots.set_vZ(i, vplus[2] + q_m * E_local[2] * 0.5 * dt);

			// Interpolate reciprocal and tangent basis vector at impurity location
			std::array<double, 9> int_rec_bas {interp_recp(bkg, slots.xidx()[i], 
				slots.yidx()[i], slots.zidx()[i], slots.x()[i], slots.y()[i], 
				slots.z()[i])};

			// Calculate velocity vector in computational coordinates using
			// interpolated reciprocal basis vector.
			slots.set_vx(i, int_rec_bas[0] * slots.vX()[i] 
				+ int_rec_bas[1] * slots.vY()[i] 
				+ int_rec_bas[2] * slots.vZ()[i]);
			slots.set_vy(i, int_rec_bas[3] * slots.vX()[i] 
				+ int_rec_bas[4] * slots.vY()[i] 
				+ int_rec_bas[5] * slots.vZ()[i]);
			slots.set_vz(i, int_rec_bas[6] * slots.vX()[i] 
				+ int_rec_bas[7] * slots.vY()[i] 
				+ int_rec_bas[8] * slots.vZ()[i]);

		}  // i loop, omp parallel for
		*/
	} // update_velocity_gpu

} // namespace Boris
