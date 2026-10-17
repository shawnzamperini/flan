/**
* @file  openadas.cu
*
* @brief CUDA routines for copying the OpenADAS rate coefficients to a GPU
*        and updating particle charge states due to ionization and
*        recombination.
*
* GPU version of the ionization/recombination routines in openadas.cpp.
* Each time step, ne and Te are interpolated at each particle's location
* (trilinear in space, linear in time), and the rate coefficients are
* bilinearly interpolated in Te and ne. Each particle can then ionize or
* recombine with a probability of rate * ne * dt.
*/
#include <cuda_runtime.h>

#include "background_device.h"
#include "openadas.h"
#include "openadas_device.h"
#include "pcg32.h"
#include "slots_device.h"


#include "interpolate.cuh"
#include "utilities.cuh"


namespace OpenADAS
{

	/**
	* Copies the OpenADAS rate coefficient data to a GPU.
	*
	* The scalar metadata (atomic number, table dimensions and charge
	* range) is copied into the returned struct, and the Te grid, ne grid
	* and flattened rate coefficient table are allocated on the device and
	* copied over. The rates are stored by charge, then Te, then ne, so
	* each charge state is a contiguous ntemp * ndens slice. The returned
	* struct must be freed with free_oa. If USE_CUDA is not defined,
	* nothing is copied and a zero-initialized struct is returned.
	*
	* @param device_id ID of the GPU to copy the data to.
	* @return Struct holding the metadata and device pointers to the data.
	*/
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


	/**
	* Frees the device memory allocated by OpenADAS::to_device.
	*
	* @param oa_d Device-side OpenADAS struct. The pointers are set to
	*        nullptr.
	* @param device_id ID of the GPU the data was copied to.
	*/
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


	/**
	* Returns a pointer to the start of the 2D (Te, ne) slice of the rate
	* coefficient table for a given charge index.
	*
	* The charge index is clamped to [0, charge_high - charge_low], so an
	* out-of-range index (e.g., -1 for recombination of a neutral) still
	* returns a valid pointer. The caller is responsible for not using the
	* rate in that case.
	*
	* @param oa_d Device-side OpenADAS struct.
	* @param charge_state Index of the charge state slice to return.
	* @return Pointer to the ntemp * ndens slice for this charge state.
	*/
	__device__ inline const double* get_charge_slice(
		const OpenADASDevice& oa_d, int charge_state)
	{
		int charge_index = charge_state;

		// Clamp to valid range
		if (charge_index < 0) charge_index = 0;
		if (charge_index > oa_d.charge_high - oa_d.charge_low)
			charge_index = oa_d.charge_high - oa_d.charge_low;

		int slice_size = oa_d.ndens * oa_d.ntemp;
		return oa_d.rates + charge_index * slice_size;
	}


	/**
	* Calculates the probabilities of a particle ionizing and recombining
	* during one time step.
	*
	* The rate coefficients are bilinearly interpolated in Te and ne, and
	* the probability is rate * ne * dt. The ionization rate for charge q
	* is at charge index q (q -> q+1), and the recombination rate is at
	* charge index q-1 (q -> q-1). A fully stripped particle (q = Z) can't
	* ionize and a neutral (q = 0) can't recombine, so those
	* probabilities are zero. The probabilities are not capped at 1.
	*
	* @param oa_ioniz_d Device-side OpenADAS struct with the ionization
	*        (SCD) rate coefficients.
	* @param oa_recomb_d Device-side OpenADAS struct with the
	*        recombination (ACD) rate coefficients.
	* @param ne Electron density at the particle location (m^-3).
	* @param te Electron temperature at the particle location (eV).
	* @param dt Time step (s).
	* @param q Charge state of the particle.
	* @param Z Atomic number of the impurity.
	* @param[out] ioniz_prob Probability of ionizing this time step.
	* @param[out] recomb_prob Probability of recombining this time step.
	*/
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
		// Te must be passed into bilinear_interp first, then ne, to match the
		// order of the indexing of the underlying 1D arrays.
		// TEMPORARY: Casting ne, Te to double just to start but should modify
		// eventually settle on float
		double ioniz_rate {can_ionize ? 
			Utilities::bilinear_interp(iz_rate_q, oa_ioniz_d.te, 
			oa_ioniz_d.ntemp, oa_ioniz_d.ne, oa_ioniz_d.ndens, 
			static_cast<double>(te), static_cast<double>(ne))
			: 0.0};
		double recomb_rate {can_recomb ? 
			Utilities::bilinear_interp(rc_rate_q, 
			oa_recomb_d.te, oa_recomb_d.ntemp, oa_recomb_d.ne, 
			oa_recomb_d.ndens,
			static_cast<double>(te), static_cast<double>(ne))
			: 0.0};

		// Compute probabilities
		float ne_dt {ne * dt};
		ioniz_prob  = ioniz_rate * ne_dt;
		recomb_prob = recomb_rate * ne_dt;
	}


	/**
	* Kernel that updates the charge state of each particle due to
	* ionization and recombination.
	*
	* One thread handles one slot, and dead particles are skipped. ne and
	* Te are interpolated at the particle location (trilinear in space
	* between the neighboring cell centers, then linear in time), and the
	* ionization and recombination probabilities are calculated with
	* calc_ioniz_recomb_probs. Each is compared to its own uniform random
	* number, so the charge changes by +1, -1 or 0 (if the particle both
	* ionizes and recombines, the two cancel). A warning counter is
	* incremented each time a probability exceeds 1, since that means dt
	* is too large.
	*
	* @param slots_d Particle slots on the device. The charge states are
	*        updated in place.
	* @param bkg_d Background plasma on the device.
	* @param oa_ioniz_d Device-side OpenADAS struct with the ionization
	*        rate coefficients.
	* @param oa_recomb_d Device-side OpenADAS struct with the
	*        recombination rate coefficients.
	* @param dt Time step (s).
	* @param rngs_d Random number generators on the device, one per slot.
	* @param d_ioniz_warnings Device counter of ionization probabilities
	*        greater than 1.
	* @param d_recomb_warnings Device counter of recombination
	*        probabilities greater than 1.
	*/
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

		// Get interpolation stencil at this location for 4D arrays
		Interpolate::InterpolationStencil stencil_4d 
			{Interpolate::build_stencil_4d(bkg_d, tidx, xidx, yidx, zidx, t, x, 
			y, z)};

		// Interpolate ne, Te at the particle location
		const double ne {Interpolate::interpolate_field_4d(bkg_d.ne,
				stencil_4d, t, x, y, z)};
		const double te {Interpolate::interpolate_field_4d(bkg_d.te,
				stencil_4d, t, x, y, z)};

		// Variables to store ionization/recombination probabilities
		float ioniz_prob {};
		float recomb_prob {};

		// Get ionization/recombination probabilities. Values are updated
		// in-place.
		calc_ioniz_recomb_probs(oa_ioniz_d, oa_recomb_d, ne, te, dt, q, Z, 
			ioniz_prob, recomb_prob);

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


	/**
	* Allocates the warning counters on a GPU and sets them to zero.
	*
	* Called once per GPU at the start of a run, not every time step. The
	* ionization and recombination counters then accumulate over every call
	* to ioniz_recomb_gpu, and the max time counter over every call to
	* ImpurityTransport::step_gpu.
	*
	* @param device_id ID of the GPU to allocate the counters on.
	* @return Struct holding the device pointers to the three counters.
	*/
	IonizRecombWarningsDevice alloc_warnings(int device_id)
	{
		IonizRecombWarningsDevice warnings_d {};
		cudaSetDevice(device_id);
		warnings_d.device_id = device_id;
		cudaMalloc(&warnings_d.max_time, sizeof(int));
		cudaMalloc(&warnings_d.ioniz, sizeof(int));
		cudaMalloc(&warnings_d.recomb, sizeof(int));
		cudaMemset(warnings_d.max_time, 0, sizeof(int));
		cudaMemset(warnings_d.ioniz, 0, sizeof(int));
		cudaMemset(warnings_d.recomb, 0, sizeof(int));
		return warnings_d;
	}


	/**
	* Copies the warning counts from the GPU and adds them to the host-side
	* counters.
	*
	* Blocks until all previously launched kernels have finished.
	*
	* @param warnings_d Device-side warning counters.
	* @param ioniz_warnings Host-side ionization warning counter.
	*        Incremented by the device count.
	* @param recomb_warnings Host-side recombination warning counter.
	*        Incremented by the device count.
	* @param max_time_ctr Host-side counter of particles that exceeded
	*        imp_max_time. Incremented by the device count.
	*/
	void add_warnings_to_host(const IonizRecombWarningsDevice& warnings_d,
		int& ioniz_warnings, int& recomb_warnings, int& max_time_ctr)
	{
		cudaSetDevice(warnings_d.device_id);

		int num_warnings {};
		cudaMemcpy(&num_warnings, warnings_d.ioniz, sizeof(int), 
			cudaMemcpyDeviceToHost);
		ioniz_warnings += num_warnings;

		num_warnings = 0;
		cudaMemcpy(&num_warnings, warnings_d.recomb, sizeof(int), 
			cudaMemcpyDeviceToHost);
		recomb_warnings += num_warnings;

		num_warnings = 0;
		cudaMemcpy(&num_warnings, warnings_d.max_time, sizeof(int), 
			cudaMemcpyDeviceToHost);
		max_time_ctr += num_warnings;
	}


	/**
	* Frees the device-side warning counters.
	*
	* @param warnings_d Device-side warning counters. The pointers are set
	*        to nullptr.
	*/
	void free_warnings(IonizRecombWarningsDevice& warnings_d)
	{
		cudaSetDevice(warnings_d.device_id);
		cudaFree(warnings_d.max_time);
		warnings_d.max_time = nullptr;
		cudaFree(warnings_d.ioniz);
		warnings_d.ioniz = nullptr;
		cudaFree(warnings_d.recomb);
		warnings_d.recomb = nullptr;
	}


	/**
	* Launches ioniz_recomb_kernel to update the charge state of every
	* particle on the GPU due to ionization and recombination.
	*
	* Uses blocks of 256 threads, one thread per slot. The warning counters
	* are allocated once with alloc_warnings and accumulate over the whole
	* run, so this is just a kernel launch with no allocations or copies.
	* The launch is asynchronous, except in DEBUG builds, where it waits
	* for the kernel to finish and prints any error. The number of
	* particles whose ionization or recombination probability exceeded 1
	* is added to the counters in warnings_d.
	*
	* @param slots_d Particle slots on the device. The charge states are
	*        updated in place.
	* @param bkg_d Background plasma on the device.
	* @param oa_ioniz_d Device-side OpenADAS struct with the ionization
	*        (SCD) rate coefficients.
	* @param oa_recomb_d Device-side OpenADAS struct with the
	*        recombination (ACD) rate coefficients.
	* @param dt Time step (s).
	* @param warnings_d Device-side warning counters from alloc_warnings.
	* @param rngs_d Random number generators on the device, one per slot.
	*/
	void ioniz_recomb_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d,
		const OpenADASDevice& oa_ioniz_d, 
		const OpenADASDevice& oa_recomb_d, const double dt,
		const IonizRecombWarningsDevice& warnings_d, pcg32* rngs_d)
	{
		// Block and grid size
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;

		ioniz_recomb_kernel<<<gridSize, blockSize>>>(slots_d, bkg_d, 
			oa_ioniz_d, oa_recomb_d, dt, rngs_d, warnings_d.ioniz, 
			warnings_d.recomb);

#ifdef DEBUG
		// Check for errors
		cudaError_t err {cudaDeviceSynchronize()};
		if (err != cudaSuccess)
			printf("ioniz_recomb_kernel error: %s\n", cudaGetErrorString(err));
#endif
	}

}  // namespace OpenADAS
