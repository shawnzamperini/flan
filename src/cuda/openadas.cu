#include <cuda_runtime.h> 

#include "background_device.h"
#include "openadas.h"
#include "openadas_device.h"
#include "pcg32.h"
#include "slots_device.h"


#include "device_constants.cuh"
#include "utilities.cuh"


namespace OpenADAS
{

	// Copy background data to device, returning an OpenADASDevice object with
	// pointers to memory location on device.
	OpenADASDevice OpenADAS::to_device(int device_id)
	{
		OpenADASDevice oa_d {};

	#ifdef USE_CUDA

		// Select device
		cudaSetDevice(device_id);
		oa_d.device_id = device_id;

		// Copy scalar metadata
		oa_d.atomic_number = m_atomic_number;
		oa_d.ndens         = m_ndens;
		oa_d.ntemp         = m_ntemp;
		oa_d.charge_low    = m_charge_low;
		oa_d.charge_high   = m_charge_high;

		const int ncharges = m_charge_high - m_charge_low + 1;
		const int nrates   = m_ndens * m_ntemp * ncharges;

		// --- Allocate device memory ---
		cudaMalloc(&oa_d.te,    m_ntemp * sizeof(double));
		cudaMalloc(&oa_d.ne,    m_ndens * sizeof(double));
		cudaMalloc(&oa_d.rates, nrates  * sizeof(double));

		// --- Copy host to device ---
		cudaMemcpy(oa_d.te, m_te.data(), m_ntemp * sizeof(double),
			cudaMemcpyHostToDevice);

		cudaMemcpy(oa_d.ne, m_ne.data(), m_ndens * sizeof(double),
			cudaMemcpyHostToDevice);

		cudaMemcpy(oa_d.rates, m_rates.get_data().data(), 
			nrates * sizeof(double), cudaMemcpyHostToDevice);

	#endif

		return oa_d;
	}  // to_device


	// Free up memory from device-side OpenADASDevice object
	void free_oa(OpenADASDevice& oa_d, int device_id)
	{
		cudaSetDevice(device_id);

		cudaFree(oa_d.te);
		oa_d.te = nullptr;

		cudaFree(oa_d.ne);
		oa_d.ne = nullptr;

		cudaFree(oa_d.rates);
		oa_d.rates = nullptr;
	}  // free_oa


	// Grab slice in OpenADASDevice.rates array the corresponds to the
	// indicated charge
	__device__ inline const double* get_charge_slice(
		const OpenADASDevice& oa_d, int charge_state)
	{
		int charge_index = charge_state - oa_d.charge_low;

		// Clamp to valid range
		if (charge_index < 0) charge_index = 0;
		if (charge_index > oa_d.charge_high - oa_d.charge_low)
			charge_index = oa_d.charge_high - oa_d.charge_low;

		int slice_size = oa_d.ndens * oa_d.ntemp;
		return oa_d.rates + charge_index * slice_size;
	}


	// Calculate the ionization/recombination probabilities at a given ne/Te
	__device__ void calc_ioniz_recomb_probs(
		const OpenADASDevice& oa_ioniz_d,
		const OpenADASDevice& oa_recomb_d, const float ne, 
		const float te, const float dt, int q, int Z, float& ioniz_prob, 
		float& recomb_prob)
	{

		// Ionization rate coefficients are indexed by charge. This is 
		// because the zeroeth charge index in the underlying rate data is
		// for neutral ionization (charge = 0), W0 --> W1+.
		const double* iz_rate_q {get_charge_slice(oa_ioniz_d, q)};

		// Recombination rate coefficients are indexed by charge-1. This is
		// because this zeroeth entry is for W1+ --> W0. So if we want that
		// rate coefficient for, say, W1+, we need to pass it charge-1 so it
		// chooses that zeroeth index.
		const double* rc_rate_q {get_charge_slice(oa_recomb_d, q-1)};

		// Branchless masks
		const bool can_ionize {q < Z};	
		const bool can_recomb {q > 0};	

		// Bilinear interpolation over the 2D slice for this charge state
		// Caller must select the correct charge slice. No point in calculating
		// if the particle can't ionize or recombine so put in an if statement.
		// TEMPORARY: Casting ne, Te to double just to start but should modify
		// eventually settle on float
		double ioniz_rate {can_ionize ? 
			Utilities::bilinear_interp(iz_rate_q, oa_ioniz_d.ne, 
			oa_ioniz_d.ndens, oa_ioniz_d.te, oa_ioniz_d.ntemp, 
			static_cast<double>(ne), static_cast<double>(te))
			: 0.0};
		double recomb_rate {can_recomb ? 
			Utilities::bilinear_interp(rc_rate_q, 
			oa_recomb_d.ne, oa_recomb_d.ndens, oa_recomb_d.te, 
			oa_recomb_d.ntemp,
			static_cast<double>(ne), static_cast<double>(te))
			: 0.0};

		// Compute probabilities
		float ne_dt {ne * dt};
		ioniz_prob  = ioniz_rate * ne_dt;
		recomb_prob = recomb_rate * ne_dt;
	}


	// Kernel for updating particle charge based on ionization/recombination
	__global__ void ioniz_recomb_kernel(Slots::SlotsDevice slots_d, 
		const Background::BackgroundDevice bkg_d, 
		OpenADASDevice oa_ioniz_d, 
		OpenADASDevice oa_recomb_d, const float dt, pcg32* rngs_d, 
		int* d_ioniz_warnings, int* d_recomb_warnings)
	{
		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;

		// Don't try and access beyond the number of slots (segfault)
		if (i >= slots_d.N) return;

		// If particle isn't alive then skip
		if (slots_d.state[i] > 0) return;

		// ---------------------------------------------------------------
		// The following algorithm is similar to that in boris.cu. It 
		// interpolates background values in time (linear) and space 
		// (trilinear). The main difference is we are just interpolating
		// different background fields (ne & Te instead of B & E). 
		// ---------------------------------------------------------------

		// Local variables
		int tidx {slots_d.tidx[i]};
		int xidx {slots_d.xidx[i]};
		int yidx {slots_d.yidx[i]};
		int zidx {slots_d.zidx[i]};

		double t {slots_d.t[i]};
		double x {slots_d.x[i]};
		double y {slots_d.y[i]};
		double z {slots_d.z[i]};

		int q {slots_d.q[i]};
		int Z {slots_d.Z};

		// Get nearest neighbor indices for each direction. These tell us
		// which direction we should interpolate towards, i.e., which
		// rectangle made by the neighboring cell centers our particle
		// is bounded by.
		// d_x is a global constant defined in cuda/device_constants.h.
		const int xidx_neighbor {Utilities::get_neighbor_index_cuda(x, 
			d_x, xidx, bkg_d.xdim)};
		const int yidx_neighbor {Utilities::get_neighbor_index_cuda(y, 
			d_y, yidx, bkg_d.ydim)};
		const int zidx_neighbor {Utilities::get_neighbor_index_cuda(z, 
			d_z, zidx, bkg_d.zdim)};

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

		// 4D -> 1D index for indexing background arrays in bkg_d. 
		// _0 is at tidx, and _1 is at tidx_neighbor
		int idx000_0_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx, xidx, yidx, zidx)};
		int idx001_0_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx, xidx, yidx, zidx_neighbor)};
		int idx010_0_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx, xidx, yidx_neighbor, zidx)};
		int idx011_0_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx, xidx, yidx_neighbor, zidx_neighbor)};
		int idx100_0_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx, xidx_neighbor, yidx, zidx)};
		int idx101_0_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx, xidx_neighbor, yidx, zidx_neighbor)};
		int idx110_0_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx, xidx_neighbor, yidx_neighbor, zidx)};
		int idx111_0_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx, xidx_neighbor, yidx_neighbor, zidx_neighbor)};

		int idx000_1_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx_neighbor, xidx, yidx, zidx)};
		int idx001_1_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx_neighbor, xidx, yidx, zidx_neighbor)};
		int idx010_1_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx_neighbor, xidx, yidx_neighbor, zidx)};
		int idx011_1_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx_neighbor, xidx, yidx_neighbor, zidx_neighbor)};
		int idx100_1_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx_neighbor, xidx_neighbor, yidx, zidx)};
		int idx101_1_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx_neighbor, xidx_neighbor, yidx, zidx_neighbor)};
		int idx110_1_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx_neighbor, xidx_neighbor, yidx_neighbor, zidx)};
		int idx111_1_4d {Utilities::calc_4d_index_cuda(bkg_d.xdim, bkg_d.ydim, 
			bkg_d.zdim, tidx_neighbor, xidx_neighbor, yidx_neighbor, 
			zidx_neighbor)};

		// Group into an array, can help reduce register pressure
		int idx0_4d[8] = {idx000_0_4d, idx001_0_4d, idx010_0_4d, idx011_0_4d,
						idx100_0_4d, idx101_0_4d, idx110_0_4d, idx111_0_4d};
		int idx1_4d[8] = {idx000_1_4d, idx001_1_4d, idx010_1_4d, idx011_1_4d,
						idx100_1_4d, idx101_1_4d, idx110_1_4d, idx111_1_4d};

		// x, y, z coordinates of two bounding vertices to interpolate 
		// between. Note these are not grid vertices, but rather are formed
		// by cell center coordinates since that's where B/E are assumed
		// to be defined. It is essentially a cell shifted by dx/2, dy/2
		// and dz/2 if that helps.
		const double t0 {d_t[tidx]};
		const double t1 {d_t[tidx_neighbor]};
		const double x0 {d_x[xidx]};
		const double x1 {d_x[xidx_neighbor]};
		const double y0 {d_y[yidx]};
		const double y1 {d_y[yidx_neighbor]};
		const double z0 {d_z[zidx]};
		const double z1 {d_z[zidx_neighbor]};

		// Distance between neighbors
		const double dx {x1 - x0};
		const double dy {y1 - y0};
		const double dz {z1 - z0};

		// Time between neighboring t frames
		const double frame_dt {t1 - t0};

		// Array to hold field component and reciprocal basis vector pointers
		// 0: ne
		// 1: Te
		// We use restrict since these arrays never overlap in memory with
		// another pointer, and allows the compiler to optimize the arrays.
		double* __restrict__ field_ptrs[2];  

		// 4D arrays
		field_ptrs[0] = bkg_d.ne;
		field_ptrs[1] = bkg_d.te;

		// Array to hold all 8 values at the vertices
		double verts0[8];
		double verts1[8];

		// This will hold the interpolated field values
		double intrp_field[2];

		// j = 0 = ne, j = 1 = Te
		for (int j = 0; j < 2; ++j)
		{

			// Calculate the 8 indices at tidx and tidx_neighbor. 
			// pragma unroll just tells the compiler to break this loop 
			// down into 8x2 individual statements. This lets it use the 
			// register more efficiently and removes loop overhead.
			#pragma unroll
			for (int l = 0; l < 8; ++l)
			{

				// 4D arrays
				verts0[l] = field_ptrs[j][idx0_4d[l]];
				verts1[l] = field_ptrs[j][idx1_4d[l]];
			}

			// Trilinear interpolation at the particle location for this
			// field component, comp0, and then comp1 is the same just
			// at the next time index. I've lined the indices of v up with
			// the corresponding ones from the trilinear_interpolate 
			// signature.
			double comp0 {Utilities::trilinear_interpolate(x0, y0, z0, 
				dx, dy, dz, verts0[0], verts0[4], verts0[2], verts0[6], 
				verts0[1], verts0[5], verts0[3], verts0[7], x, y, z)};

			double comp1 {Utilities::trilinear_interpolate(x0, y0, z0, 
				dx, dy, dz, verts1[0], verts1[4], verts1[2], verts1[6], 
				verts1[1], verts1[5], verts1[3], verts1[7], x, y, z)};

			// Now linearly interpolate in time (this is just point slope)
			const double slope {(comp1 - comp0) / frame_dt};
			intrp_field[j] = slope * (t - t0) + comp0;
		}

		// Variables to store ionization/recombination probabilities
		float ioniz_prob {};
		float recomb_prob {};

		// Get ionization/recombination probabilities. Values are updated
		// in-place. intrp_field[0] = ne, intrp_field[1] = Te.
		calc_ioniz_recomb_probs(oa_ioniz_d, oa_recomb_d, intrp_field[0], 
			intrp_field[1], dt, q, Z, ioniz_prob, recomb_prob);

		// Increment warnings if either probability is greater than 1. We
		// could easily do this without if statements, but in this case these
		// warnings should be rare (in any valid simulation at least) so this
		// method avoids the atomicAdds and only executes them when needed.
		if (ioniz_prob > 1.0)
			atomicAdd(d_ioniz_warnings, 1);
		if (recomb_prob > 1.0)
			atomicAdd(d_recomb_warnings, 1);

		// Pull random numbers
		const double r1 = rngs_d[i].next_double();
		const double r2 = rngs_d[i].next_double();

		// Compare to probs and update charge state if
		const int ionize = (r1 < ioniz_prob);
		const int recomb = (r2 < recomb_prob);

		// Clever branchless way to adjust the charge
		const int dq = ionize - recomb;
		slots_d.q[i] = q + dq;
	}


	// Wrapper to call ionization/recombination kernel
	void ioniz_recomb_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d,
		const OpenADASDevice& oa_ioniz_d, 
		const OpenADASDevice& oa_recomb_d, const double dt,
		int& ioniz_warnings, int& recomb_warnings, pcg32* rngs_d)
	{

		// Counter for ionization/recombination warnings
		int* d_ioniz_warnings;
		int* d_recomb_warnings;
		cudaMalloc(&d_ioniz_warnings, sizeof(int));
		cudaMalloc(&d_recomb_warnings, sizeof(int));
		cudaMemset(d_ioniz_warnings, 0, sizeof(int));
		cudaMemset(d_recomb_warnings, 0, sizeof(int));

		// Block and grid size
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;

		ioniz_recomb_kernel<<<gridSize, blockSize>>>(slots_d, bkg_d, 
			oa_ioniz_d, oa_recomb_d, dt, rngs_d, d_ioniz_warnings, 
			d_recomb_warnings);

		// Retrieve how many times the warning was incremented for each
		int num_warnings = 0;
		cudaMemcpy(&num_warnings, d_ioniz_warnings, sizeof(int), 
			cudaMemcpyDeviceToHost);
		ioniz_warnings += num_warnings;

		num_warnings = 0;
		cudaMemcpy(&num_warnings, d_recomb_warnings, sizeof(int), 
			cudaMemcpyDeviceToHost);
		recomb_warnings += num_warnings;
	}


}  // namespace OpenADAS
