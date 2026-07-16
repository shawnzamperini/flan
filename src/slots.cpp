#include <iostream>
#include <vector>

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

	// Copy data to device and return SlotsDevice struct
	SlotsDevice Slots::to_device(int device_id)
	{

		// Create struct to hold device-side memory locations
		SlotsDevice slots_d {};

#ifdef USE_CUDA

		slots_d.device_id = device_id;
		slots_d.N = m_N;

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

#else
		std::cerr << "Error! Slots.to_device() was called but GPU support"
			<< " was not compiled in.\n";
#endif

		return slots_d;
	}

	// Free up memory from device-side SlotsDevice object
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
#endif
		return;
	}

	// Initialize a new particle and return it
	ParticleInit make_new_particle() 
	{
		ParticleInit p;

		// To-do: Replace with actual initialization logic
		p.t = 0.0;
		p.x = 0.0;
		p.y = 0.0;
		p.z = 0.0;

		p.tidx = 0;
		p.xidx = 0;
		p.yidx = 0;
		p.zidx = 0;

		p.vx = 0.0;
		p.vy = 0.0;
		p.vz = 0.0;

		p.vX = 0.0;
		p.vY = 0.0;
		p.vZ = 0.0;
		
		p.weight = 1.0;
		p.q = 1;

		return p;
	}

	// Replace dead particles with alive ones, as long as remaining particles
	// are greater than zero.
	void fill_slots_cpu(Slots& slots, int& rem_parts)
	{
		// Loop through each slot as long as there are remaining particles
		// Make this an OpenMP loop, watch for race condition on rem_parts
		for (int i = 0; i < slots.N() && rem_parts > 0; i++)
		{
			// If particle is dead
			if (slots.state()[i] > 0)
			{
				ParticleInit p = make_new_particle();

				// Assign initial positions and velocities
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

				// Reassign to alive, decrease remaining particles by one
				slots.set_state(i, 0);
				rem_parts--;
			}
		}
	}

}  // namespace Slots
