#include <iostream>
#include <omp.h>
#include <vector>

#include "background.h"
#include "pcg32.h"
#include "slots.h"
#include "slots_device.h"

#ifdef USE_CUDA
#include <cuda_runtime.h> 
#endif

namespace Slots
{
	// Constructor
	Slots::Slots(int N)
		: m_N {N}
	{
		// Resize particle vectors
		m_t.resize(N);
		m_x.resize(N);
		m_y.resize(N);
		m_z.resize(N);
		m_vx.resize(N);
		m_vy.resize(N);
		m_vz.resize(N);
		m_vX.resize(N);
		m_vY.resize(N);
		m_vZ.resize(N);
		m_tidx.resize(N);
		m_xidx.resize(N);
		m_yidx.resize(N);
		m_zidx.resize(N);
		m_weight.resize(N);
		m_q.resize(N);
		m_state.resize(N);
	}

	// Destructor not actually needed because we use std::vector, which
	// free itself automatically. GPU memory must be freed with free_slots.
	Slots::~Slots() = default;
	
	// Accessors
	int Slots::N() const noexcept {return m_N;}
	int Slots::Z() const noexcept {return m_Z;}
	double Slots::mass() const noexcept {return m_mass;}
	const std::vector<double>& Slots::t() const noexcept {return m_t;}
	const std::vector<double>& Slots::x() const noexcept {return m_x;}
	const std::vector<double>& Slots::y() const noexcept {return m_y;}
	const std::vector<double>& Slots::z() const noexcept {return m_z;}
	const std::vector<double>& Slots::vx() const noexcept {return m_vx;}
	const std::vector<double>& Slots::vy() const noexcept {return m_vy;}
	const std::vector<double>& Slots::vz() const noexcept {return m_vz;}
	const std::vector<double>& Slots::vX() const noexcept {return m_vX;}
	const std::vector<double>& Slots::vY() const noexcept {return m_vY;}
	const std::vector<double>& Slots::vZ() const noexcept {return m_vZ;}
	const std::vector<int>& Slots::tidx() const noexcept {return m_tidx;}
	const std::vector<int>& Slots::xidx() const noexcept {return m_xidx;}
	const std::vector<int>& Slots::yidx() const noexcept {return m_yidx;}
	const std::vector<int>& Slots::zidx() const noexcept {return m_zidx;}
	const std::vector<double>& Slots::weight() const noexcept {return m_weight;}
	const std::vector<int>& Slots::q() const noexcept {return m_q;}
	const std::vector<int>& Slots::state() const noexcept {return m_state;}

	void Slots::set_mass(double mass) {m_mass = mass;}
	void Slots::set_Z(int Z) {m_Z = Z;}

	// Element level setters
	void Slots::set_t(int i, double val) {m_t[i] = val;}
	void Slots::set_x(int i, double val) {m_x[i] = val;}
	void Slots::set_y(int i, double val) {m_y[i] = val;}
	void Slots::set_z(int i, double val) {m_z[i] = val;}
	void Slots::set_vx(int i, double val) {m_vx[i] = val;}
	void Slots::set_vy(int i, double val) {m_vy[i] = val;}
	void Slots::set_vz(int i, double val) {m_vz[i] = val;}
	void Slots::set_vX(int i, double val) {m_vX[i] = val;}
	void Slots::set_vY(int i, double val) {m_vY[i] = val;}
	void Slots::set_vZ(int i, double val) {m_vZ[i] = val;}
	void Slots::set_tidx(int i, int val) {m_tidx[i] = val;}
	void Slots::set_xidx(int i, int val) {m_xidx[i] = val;}
	void Slots::set_yidx(int i, int val) {m_yidx[i] = val;}
	void Slots::set_zidx(int i, int val) {m_zidx[i] = val;}
	void Slots::set_weight(int i, double val) {m_weight[i] = val;}
	void Slots::set_q(int i, int val) {m_q[i] = val;}
	void Slots::set_state(int i, int val) {m_state[i] = val;}

	// I think this should probably be moved to slots.cu...
	// Copy data to device and return SlotsDevice struct
	SlotsDevice Slots::to_device(int device_id)
	{

		// Create struct on the host that holds device-side memory locations.
		SlotsDevice slots_d {};

#ifdef USE_CUDA

		// These are just scalars, no allocation needed
		slots_d.device_id = device_id;
		slots_d.N = m_N;
		slots_d.Z = m_Z;
		slots_d.mass = m_mass;

		// GPU allocation
		cudaSetDevice(device_id);

		cudaMalloc(&slots_d.t, m_N * sizeof(double));
		cudaMalloc(&slots_d.x, m_N * sizeof(double));
		cudaMalloc(&slots_d.y, m_N * sizeof(double));
		cudaMalloc(&slots_d.z, m_N * sizeof(double));
		cudaMalloc(&slots_d.tidx, m_N * sizeof(int));
		cudaMalloc(&slots_d.xidx, m_N * sizeof(int));
		cudaMalloc(&slots_d.yidx, m_N * sizeof(int));
		cudaMalloc(&slots_d.zidx, m_N * sizeof(int));
		cudaMalloc(&slots_d.vx, m_N * sizeof(double));
		cudaMalloc(&slots_d.vy, m_N * sizeof(double));
		cudaMalloc(&slots_d.vz, m_N * sizeof(double));
		cudaMalloc(&slots_d.vX, m_N * sizeof(double));
		cudaMalloc(&slots_d.vY, m_N * sizeof(double));
		cudaMalloc(&slots_d.vZ, m_N * sizeof(double));
		cudaMalloc(&slots_d.weight, m_N * sizeof(double));
		cudaMalloc(&slots_d.q, m_N * sizeof(int));
		cudaMalloc(&slots_d.state, m_N * sizeof(int));
		cudaMalloc(&slots_d.all_dead, sizeof(bool));
		cudaMalloc(&slots_d.counter, sizeof(int));
		cudaMalloc(&slots_d.alive, sizeof(int));
		cudaMalloc(&slots_d.num_dead, sizeof(int));
		cudaMalloc(&slots_d.dead_indices, m_N * sizeof(int));

		// Copy to device
		cudaMemcpy(slots_d.t, m_t.data(), m_N * sizeof(double), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(slots_d.x, m_x.data(), m_N * sizeof(double), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(slots_d.y, m_y.data(), m_N * sizeof(double), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(slots_d.z, m_z.data(), m_N * sizeof(double), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(slots_d.tidx, m_tidx.data(), m_N * sizeof(int), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(slots_d.xidx, m_xidx.data(), m_N * sizeof(int), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(slots_d.yidx, m_yidx.data(), m_N * sizeof(int), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(slots_d.zidx, m_zidx.data(), m_N * sizeof(int), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(slots_d.vx, m_vx.data(), m_N * sizeof(double),
			cudaMemcpyHostToDevice);
		cudaMemcpy(slots_d.vy, m_vy.data(), m_N * sizeof(double),
			cudaMemcpyHostToDevice);
		cudaMemcpy(slots_d.vz, m_vz.data(), m_N * sizeof(double),
			cudaMemcpyHostToDevice);
		cudaMemcpy(slots_d.vX, m_vX.data(), m_N * sizeof(double),
			cudaMemcpyHostToDevice);
		cudaMemcpy(slots_d.vY, m_vY.data(), m_N * sizeof(double),
			cudaMemcpyHostToDevice);
		cudaMemcpy(slots_d.vZ, m_vZ.data(), m_N * sizeof(double),
			cudaMemcpyHostToDevice);
		cudaMemcpy(slots_d.weight, m_weight.data(), m_N * sizeof(double),
			cudaMemcpyHostToDevice);
		cudaMemcpy(slots_d.q,  m_q.data(), m_N * sizeof(int),
			cudaMemcpyHostToDevice);
		cudaMemcpy(slots_d.state,  m_state.data(), m_N * sizeof(int),
			cudaMemcpyHostToDevice);

		// Initialize all_dead to false
		bool init = false;
		cudaMemcpy(slots_d.all_dead, &init, sizeof(bool), 
			cudaMemcpyHostToDevice);

#else
		std::cerr << "Error! Slots.to_device() was called but GPU support"
			<< " was not compiled in.\n";
#endif

		return slots_d;
	}

	// Free up memory from device-side SlotsDevice object
	Slots Slots::to_host(const SlotsDevice& slots_d)
	{
		Slots slots(slots_d.N);
		slots.set_mass(slots_d.mass);

#ifdef USE_CUDA
		cudaSetDevice(slots_d.device_id);
		cudaMemcpy(slots.m_t.data(), slots_d.t, slots_d.N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		cudaMemcpy(slots.m_x.data(), slots_d.x, slots_d.N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		cudaMemcpy(slots.m_y.data(), slots_d.y, slots_d.N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		cudaMemcpy(slots.m_z.data(), slots_d.z, slots_d.N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		cudaMemcpy(slots.m_vx.data(), slots_d.vx, slots_d.N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		cudaMemcpy(slots.m_vy.data(), slots_d.vy, slots_d.N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		cudaMemcpy(slots.m_vz.data(), slots_d.vz, slots_d.N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		cudaMemcpy(slots.m_vX.data(), slots_d.vX, slots_d.N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		cudaMemcpy(slots.m_vY.data(), slots_d.vY, slots_d.N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		cudaMemcpy(slots.m_vZ.data(), slots_d.vZ, slots_d.N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		cudaMemcpy(slots.m_tidx.data(), slots_d.tidx, slots_d.N * sizeof(int),	
			cudaMemcpyDeviceToHost);
		cudaMemcpy(slots.m_xidx.data(), slots_d.xidx, slots_d.N * sizeof(int), 
			cudaMemcpyDeviceToHost);
		cudaMemcpy(slots.m_yidx.data(), slots_d.yidx, slots_d.N * sizeof(int), 
			cudaMemcpyDeviceToHost);
		cudaMemcpy(slots.m_zidx.data(), slots_d.zidx, slots_d.N * sizeof(int), 
			cudaMemcpyDeviceToHost);
		cudaMemcpy(slots.m_weight.data(), slots_d.weight, 
			slots_d.N * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(slots.m_q.data(), slots_d.q, slots_d.N * sizeof(int), 
			cudaMemcpyDeviceToHost);
		cudaMemcpy(slots.m_state.data(), slots_d.state, 
			slots_d.N * sizeof(int), cudaMemcpyDeviceToHost);
#endif
		return slots;
	}

	void free_slots(SlotsDevice& slots_d, int device_id)
	{

#ifdef USE_CUDA

		// Free up memory on GPU
		cudaSetDevice(device_id);

		cudaFree(slots_d.t);
		cudaFree(slots_d.x);
		cudaFree(slots_d.y);
		cudaFree(slots_d.z);
		cudaFree(slots_d.tidx);
		cudaFree(slots_d.xidx);
		cudaFree(slots_d.yidx);
		cudaFree(slots_d.zidx);
		cudaFree(slots_d.vx);
		cudaFree(slots_d.vy);
		cudaFree(slots_d.vz);
		cudaFree(slots_d.vX);
		cudaFree(slots_d.vY);
		cudaFree(slots_d.vZ);
		cudaFree(slots_d.weight);
		cudaFree(slots_d.q);
		cudaFree(slots_d.state);
		cudaFree(slots_d.all_dead);
		cudaFree(slots_d.counter);
		cudaFree(slots_d.alive);
		cudaFree(slots_d.num_dead);
		cudaFree(slots_d.dead_indices);

		slots_d.t = nullptr;
		slots_d.x = nullptr;
		slots_d.y = nullptr;
		slots_d.z = nullptr;
		slots_d.tidx = nullptr;
		slots_d.xidx = nullptr;
		slots_d.yidx = nullptr;
		slots_d.zidx = nullptr;
		slots_d.vx = nullptr;
		slots_d.vy = nullptr;
		slots_d.vz = nullptr;
		slots_d.vX = nullptr;
		slots_d.vY = nullptr;
		slots_d.vZ = nullptr;
		slots_d.weight = nullptr;
		slots_d.q = nullptr;
		slots_d.state = nullptr;
		slots_d.all_dead = nullptr;
#endif
		return;
	}

	// Function to decide starting t,x,y,z based on input options
	double get_birth_val(const int start_opt_int, const double start_val, 
		const double range_min, const double range_max, const double bkg_min, 
		const double bkg_max, pcg32& rng)
	{
		// Start at specific point
		double return_start_val {};
		if (start_opt_int == 0)
		{
			 return_start_val = start_val;
		}

		// Start between a user-specified range 
		else if (start_opt_int == 1)
		{
			return_start_val = range_min + (range_max - range_min) 
				* rng.next_double();
		}

		// Start between the full range of the simulation volume 
		else if (start_opt_int == 2)
		{
			return_start_val = bkg_min + (bkg_max - bkg_min) 
				* rng.next_double();
		}
		
		return return_start_val;
	}

	// Initialize a new particle and return it
	ParticleInit make_new_particle(const Background::Background& bkg,
		const Options::Options& opts, pcg32& rng) 
	{
		ParticleInit p;

		// Initialize starting time/location based on input options
		p.t = get_birth_val(opts.imp_tstart_opt_int(), 
			opts.imp_tstart_val(), opts.imp_trange_min(), opts.imp_trange_max(), 
			bkg.get_t_min(), bkg.get_t_max(), rng);

		p.x = get_birth_val(opts.imp_xstart_opt_int(), 
			opts.imp_xstart_val(), opts.imp_xrange_min(), opts.imp_xrange_max(), 
			bkg.get_x_min(), bkg.get_x_max(), rng);

		p.y = get_birth_val(opts.imp_ystart_opt_int(), 
			opts.imp_ystart_val(), opts.imp_yrange_min(), opts.imp_yrange_max(), 
			bkg.get_y_min(), bkg.get_y_max(), rng);

		p.z = get_birth_val(opts.imp_zstart_opt_int(), 
			opts.imp_zstart_val(), opts.imp_zrange_min(), opts.imp_zrange_max(), 
			bkg.get_z_min(), bkg.get_z_max(), rng);


		// These get assigned in the main_loop right after fill_slots
		p.tidx = 0;
		p.xidx = 0;
		p.yidx = 0;
		p.zidx = 0;

		// Still need to implement the normal logic here
		p.vx = 0.0;
		p.vy = 0.0;
		p.vz = 0.0;

		p.vX = 0.0;
		p.vY = 5000.0;
		p.vZ = 0.0;
		
		p.weight = 1.0;
		p.q = opts.imp_init_charge(); 

		return p;
	}

	// Replace dead particles with alive ones, as long as remaining particles
	// are greater than zero.
	void fill_slots_cpu(Slots& slots, int& rem_parts, int& alive_slots,
		const Background::Background& bkg, const Options::Options& opts,
		std::vector<pcg32>& rngs)
	{
		int N = slots.N();
		int dead_count = 0;

		// First pass: count dead slots (reduction)
		#pragma omp parallel for reduction(+:dead_count)
		for (int i = 0; i < N; i++)
		{
			if (slots.state()[i] > 0) dead_count++;
		}

		// Determine how many particles to revive, capping off at rem_parts
		// if there are more dead particles than remaining.
		int revive = (rem_parts < dead_count) ? rem_parts : dead_count;
		if (revive < 0) revive = 0;
		rem_parts -= revive;

		// Second pass: revive without atomics
		int revived = 0;
		#pragma omp parallel
		{

			// Grab our RNG for this thread, each thread has it own (and they're
			// seeded uniquely).
			int tid = omp_get_thread_num();
			pcg32& rng = rngs[tid];

			#pragma omp for
			for (int i = 0; i < N; i++)
			{
				if (slots.state()[i] > 0)
				{
					int claim = 0;
					#pragma omp atomic capture
					{
						claim = revived;
						revived++;
					}
					if (claim >= revive)
						continue;

					ParticleInit p = make_new_particle(bkg, opts, rng);
					slots.set_t(i, p.t);
					slots.set_x(i, p.x);
					slots.set_y(i, p.y);
					slots.set_z(i, p.z);
					slots.set_tidx(i, p.tidx);
					slots.set_xidx(i, p.xidx);
					slots.set_yidx(i, p.yidx);
					slots.set_zidx(i, p.zidx);
					slots.set_vx(i, p.vx);
					slots.set_vy(i, p.vy);
					slots.set_vz(i, p.vz);
					slots.set_vX(i, p.vX);
					slots.set_vY(i, p.vY);
					slots.set_vZ(i, p.vZ);
					slots.set_q(i, p.q);
					slots.set_weight(i, 1.0);
					slots.set_state(i, 0);
				}
			}  // i loop
		}  // pragma omp parallel

		if (rem_parts < 0) rem_parts = 0;

		// Number of alive slots is the number of zeros
		alive_slots = std::count(slots.state().begin(), slots.state().end(), 0);
	}

}  // namespace Slots
