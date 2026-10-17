/**
* @file  slots.cpp
*
* @brief Storage for the particles being followed and routines to create new
*        particles and move them to and from the GPU.
*
* A Slots object holds a fixed number of particle "slots". Each slot holds one
* particle that is followed until it dies, at which point the slot is refilled
* with a new particle by fill_slots_cpu (or its CUDA version) until the
* requested number of particles has been followed. The reason for this layout
* is to enable a structure of arrays (SoA) memeory layout, which is more
* efficient memory-wise.
*/
#include <algorithm>
#include <array>
#include <iostream>
#include <omp.h>
#include <vector>

#include "background.h"
#include "indices.h"
#include "pcg32.h"
#include "slots.h"
#include "slots_device.h"

#ifdef USE_CUDA
#include <cuda_runtime.h> 
#endif

namespace Slots
{
	/**
	* Creates a Slots object with room for N particles.
	*
	* Every per-particle vector is resized to N. The values are zero
	* initialized, so every slot starts with state 0 (alive). The mass and
	* atomic number are not set here, use set_mass and set_Z.
	*
	* @param N Number of particle slots.
	*/
	Slots::Slots(int N)
		: m_N {N}
	{
		// Resize particle vectors
		m_pt.resize(N);
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
		m_iters.resize(N);
		m_tidx.resize(N);
		m_xidx.resize(N);
		m_yidx.resize(N);
		m_zidx.resize(N);
		m_weight.resize(N);
		m_q.resize(N);
		m_state.resize(N);
	}

	/**
	* Destroys the Slots object.
	*
	* Not actually needed because we use std::vector, which frees itself
	* automatically. GPU memory must be freed with free_slots.
	*/
	Slots::~Slots() = default;
	
	/**
	* Returns the number of particle slots.
	*
	* @return Number of slots.
	*/
	int Slots::N() const noexcept {return m_N;}

	/**
	* Returns the atomic number of the impurity. All particles share the
	* same atomic number.
	*
	* @return Atomic number.
	*/
	int Slots::Z() const noexcept {return m_Z;}

	/**
	* Returns the mass of the impurity. All particles share the same mass.
	*
	* @return Impurity mass (kg).
	*/
	double Slots::mass() const noexcept {return m_mass;}

	/**
	* Returns the time each particle has been followed since it was born (s).
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<double>& Slots::pt() const noexcept {return m_pt;}

	/**
	* Returns the simulation time of each particle (s).
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<double>& Slots::t() const noexcept {return m_t;}

	/**
	* Returns the computational x coordinate of each particle.
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<double>& Slots::x() const noexcept {return m_x;}

	/**
	* Returns the computational y coordinate of each particle.
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<double>& Slots::y() const noexcept {return m_y;}

	/**
	* Returns the computational z coordinate of each particle.
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<double>& Slots::z() const noexcept {return m_z;}

	/**
	* Returns the velocity of each particle along computational x.
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<double>& Slots::vx() const noexcept {return m_vx;}

	/**
	* Returns the velocity of each particle along computational y.
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<double>& Slots::vy() const noexcept {return m_vy;}

	/**
	* Returns the velocity of each particle along computational z.
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<double>& Slots::vz() const noexcept {return m_vz;}

	/**
	* Returns the Cartesian X velocity of each particle (m/s).
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<double>& Slots::vX() const noexcept {return m_vX;}

	/**
	* Returns the Cartesian Y velocity of each particle (m/s).
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<double>& Slots::vY() const noexcept {return m_vY;}

	/**
	* Returns the Cartesian Z velocity of each particle (m/s).
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<double>& Slots::vZ() const noexcept {return m_vZ;}

	/**
	* Returns the number of time steps each particle has taken.
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<int>& Slots::iters() const noexcept {return m_iters;}

	/**
	* Returns the index of the background time frame nearest to each particle.
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<int>& Slots::tidx() const noexcept {return m_tidx;}

	/**
	* Returns the x index of the grid cell containing each particle.
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<int>& Slots::xidx() const noexcept {return m_xidx;}

	/**
	* Returns the y index of the grid cell containing each particle.
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<int>& Slots::yidx() const noexcept {return m_yidx;}

	/**
	* Returns the z index of the grid cell containing each particle.
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<int>& Slots::zidx() const noexcept {return m_zidx;}

	/**
	* Returns the Monte Carlo weight of each particle.
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<double>& Slots::weight() const noexcept {return m_weight;}

	/**
	* Returns the charge state of each particle.
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<int>& Slots::q() const noexcept {return m_q;}

	/**
	* Returns the state of each particle, where 0 is alive and anything
	* greater than 0 is dead.
	*
	* @return Vector of length N() with one entry per slot.
	*/
	const std::vector<int>& Slots::state() const noexcept {return m_state;}

	/**
	* Sets the mass of the impurity, shared by all particles.
	*
	* @param mass Impurity mass (kg).
	*/
	void Slots::set_mass(double mass) {m_mass = mass;}

	/**
	* Sets the atomic number of the impurity, shared by all particles.
	*
	* @param Z Atomic number.
	*/
	void Slots::set_Z(int Z) {m_Z = Z;}

	/**
	* Sets the pt value of a single slot.
	*
	* @param i Index of the slot.
	* @param val Time since birth to store in slot i.
	*/
	void Slots::set_pt(int i, double val) {m_pt[i] = val;}

	/**
	* Sets the t value of a single slot.
	*
	* @param i Index of the slot.
	* @param val Simulation time to store in slot i.
	*/
	void Slots::set_t(int i, double val) {m_t[i] = val;}

	/**
	* Sets the x value of a single slot.
	*
	* @param i Index of the slot.
	* @param val Computational x coordinate to store in slot i.
	*/
	void Slots::set_x(int i, double val) {m_x[i] = val;}

	/**
	* Sets the y value of a single slot.
	*
	* @param i Index of the slot.
	* @param val Computational y coordinate to store in slot i.
	*/
	void Slots::set_y(int i, double val) {m_y[i] = val;}

	/**
	* Sets the z value of a single slot.
	*
	* @param i Index of the slot.
	* @param val Computational z coordinate to store in slot i.
	*/
	void Slots::set_z(int i, double val) {m_z[i] = val;}

	/**
	* Sets the vx value of a single slot.
	*
	* @param i Index of the slot.
	* @param val Velocity along computational x to store in slot i.
	*/
	void Slots::set_vx(int i, double val) {m_vx[i] = val;}

	/**
	* Sets the vy value of a single slot.
	*
	* @param i Index of the slot.
	* @param val Velocity along computational y to store in slot i.
	*/
	void Slots::set_vy(int i, double val) {m_vy[i] = val;}

	/**
	* Sets the vz value of a single slot.
	*
	* @param i Index of the slot.
	* @param val Velocity along computational z to store in slot i.
	*/
	void Slots::set_vz(int i, double val) {m_vz[i] = val;}

	/**
	* Sets the vX value of a single slot.
	*
	* @param i Index of the slot.
	* @param val Cartesian X velocity (m/s) to store in slot i.
	*/
	void Slots::set_vX(int i, double val) {m_vX[i] = val;}

	/**
	* Sets the vY value of a single slot.
	*
	* @param i Index of the slot.
	* @param val Cartesian Y velocity (m/s) to store in slot i.
	*/
	void Slots::set_vY(int i, double val) {m_vY[i] = val;}

	/**
	* Sets the vZ value of a single slot.
	*
	* @param i Index of the slot.
	* @param val Cartesian Z velocity (m/s) to store in slot i.
	*/
	void Slots::set_vZ(int i, double val) {m_vZ[i] = val;}

	/**
	* Sets the iters value of a single slot.
	*
	* @param i Index of the slot.
	* @param val Number of time steps taken to store in slot i.
	*/
	void Slots::set_iters(int i, int val) {m_iters[i] = val;}

	/**
	* Sets the tidx value of a single slot.
	*
	* @param i Index of the slot.
	* @param val Index of the nearest background time frame to store in slot i.
	*/
	void Slots::set_tidx(int i, int val) {m_tidx[i] = val;}

	/**
	* Sets the xidx value of a single slot.
	*
	* @param i Index of the slot.
	* @param val x index of the containing grid cell to store in slot i.
	*/
	void Slots::set_xidx(int i, int val) {m_xidx[i] = val;}

	/**
	* Sets the yidx value of a single slot.
	*
	* @param i Index of the slot.
	* @param val y index of the containing grid cell to store in slot i.
	*/
	void Slots::set_yidx(int i, int val) {m_yidx[i] = val;}

	/**
	* Sets the zidx value of a single slot.
	*
	* @param i Index of the slot.
	* @param val z index of the containing grid cell to store in slot i.
	*/
	void Slots::set_zidx(int i, int val) {m_zidx[i] = val;}

	/**
	* Sets the weight value of a single slot.
	*
	* @param i Index of the slot.
	* @param val Monte Carlo weight to store in slot i.
	*/
	void Slots::set_weight(int i, double val) {m_weight[i] = val;}

	/**
	* Sets the q value of a single slot.
	*
	* @param i Index of the slot.
	* @param val Charge state to store in slot i.
	*/
	void Slots::set_q(int i, int val) {m_q[i] = val;}

	/**
	* Sets the state value of a single slot.
	*
	* @param i Index of the slot.
	* @param val State to store in slot i, where 0 is alive and anything
	*        greater than 0 is dead.
	*/
	void Slots::set_state(int i, int val) {m_state[i] = val;}

	/**
	* Allocates GPU memory for the particles and copies them to the device.
	*
	* I think this should probably be moved to slots.cu...
	*
	* The per-particle arrays are allocated and copied, except iters,
	* which are not stored on the device. The work arrays used by the GPU
	* kernels (all_dead, num_dead and dead_indices) are also
	* allocated, and all_dead is initialized to false. If GPU support was not
	* compiled in, an error is printed and an empty SlotsDevice is returned.
	*
	* @param device_id ID of the GPU to allocate the memory on.
	* @return SlotsDevice struct holding the device pointers. The memory must
	*         be released with free_slots.
	*/
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

		cudaMalloc(&slots_d.pt, m_N * sizeof(double));
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
		cudaMalloc(&slots_d.iters, m_N * sizeof(int));
		cudaMalloc(&slots_d.weight, m_N * sizeof(double));
		cudaMalloc(&slots_d.q, m_N * sizeof(int));
		cudaMalloc(&slots_d.state, m_N * sizeof(int));
		cudaMalloc(&slots_d.all_dead, sizeof(bool));
		cudaMalloc(&slots_d.num_dead, sizeof(int));
		cudaMalloc(&slots_d.dead_indices, m_N * sizeof(int));

		// Copy to device
		cudaMemcpy(slots_d.pt, m_pt.data(), m_N * sizeof(double), 
			cudaMemcpyHostToDevice);
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

	/**
	* Copies the particles from the device back into a new host Slots object.
	*
	* The mass and the per-particle arrays on the device are copied back.
	* iters and the atomic number are not stored on the device, so they are
	* left at zero in the returned object. The device memory is not freed,
	* use free_slots for that. If GPU support was not compiled in, only N and
	* the mass are set.
	*
	* @param slots_d SlotsDevice struct returned by to_device.
	* @return Host Slots object holding a copy of the device data.
	*/
	Slots Slots::to_host(const SlotsDevice& slots_d)
	{
		Slots slots(slots_d.N);
		slots.set_mass(slots_d.mass);

#ifdef USE_CUDA
		cudaSetDevice(slots_d.device_id);
		cudaMemcpy(slots.m_pt.data(), slots_d.pt, slots_d.N * sizeof(double), 
			cudaMemcpyDeviceToHost);
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
		cudaMemcpy(slots.m_iters.data(), slots_d.iters, slots_d.N * sizeof(int),	
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

	/**
	* Frees the GPU memory held by a SlotsDevice struct.
	*
	* Every device array allocated in Slots::to_device is freed. The particle
	* array pointers and all_dead are then set to nullptr. Does nothing if GPU
	* support was not compiled in.
	*
	* @param slots_d SlotsDevice struct whose memory is freed.
	* @param device_id ID of the GPU the memory was allocated on.
	*/
	void free_slots(SlotsDevice& slots_d, int device_id)
	{

#ifdef USE_CUDA

		// Free up memory on GPU
		cudaSetDevice(device_id);

		cudaFree(slots_d.pt);
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
		cudaFree(slots_d.iters);
		cudaFree(slots_d.weight);
		cudaFree(slots_d.q);
		cudaFree(slots_d.state);
		cudaFree(slots_d.all_dead);
		cudaFree(slots_d.num_dead);
		cudaFree(slots_d.dead_indices);

		slots_d.pt = nullptr;
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
		slots_d.iters = nullptr;
		slots_d.weight = nullptr;
		slots_d.q = nullptr;
		slots_d.state = nullptr;
		slots_d.all_dead = nullptr;
		slots_d.num_dead = nullptr;
		slots_d.dead_indices = nullptr;
#endif
		return;
	}

	/**
	* Picks the starting value of t, x, y or z for a new particle based on
	* the input options.
	*
	* The same function is used for each coordinate, with the matching
	* imp_*start_opt option deciding where the particle starts:
	*   0: At start_val.
	*   1: Uniformly random between range_min and range_max.
	*   2: Uniformly random over the full simulation volume, between bkg_min
	*      and bkg_max.
	* Any other value returns 0.
	*
	* @param start_opt_int Control integer of the imp_*start_opt option.
	* @param start_val Starting value used when start_opt_int is 0.
	* @param range_min Lower end of the range used when start_opt_int is 1.
	* @param range_max Upper end of the range used when start_opt_int is 1.
	* @param bkg_min Minimum of the simulation volume, used when start_opt_int
	*        is 2.
	* @param bkg_max Maximum of the simulation volume, used when start_opt_int
	*        is 2.
	* @param rng Random number generator for this thread.
	* @return Starting value of the coordinate.
	*/
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


	/**
	* Picks the starting Cartesian velocity of a new particle.
	*
	* For the test backgrounds (bkg_source = "test") the velocity is fixed
	* by the test case: 10,000 m/s in Y for the gyrate, ExB, grad-B and
	* polarization tests, 30,000 m/s along the field line (in the X-Y plane,
	* at angle z) for the curvature test, and at rest for the friction test.
	*
	* Otherwise each component is sampled from a Maxwellian with zero mean
	* and sigma = sqrt(kT/m), giving an isotropic velocity distribution. The
	* temperature is imp_temp_start_val if imp_temp_start_opt is
	* "single_value", or the main ion temperature interpolated at the
	* particle location with Background::interp_ti if it is "main_ion".
	* The interpolated ion temperature is held at or above min_ti, since it
	* can go negative near the grid edges.
	*
	* @param bkg Background plasma, used for the ion temperature.
	* @param opts Simulation options.
	* @param t Starting time of the particle (s).
	* @param x Starting computational x coordinate.
	* @param y Starting computational y coordinate.
	* @param z Starting computational z coordinate. In the curvature test,
	*        z is the toroidal angle.
	* @param tidx Index of the nearest background time frame. Currently
	*        unused.
	* @param xidx x index of the grid cell containing the particle.
	*        Currently unused.
	* @param yidx y index of the grid cell containing the particle.
	*        Currently unused.
	* @param zidx z index of the grid cell containing the particle.
	*        Currently unused.
	* @param rng Random number generator for this thread.
	* @return Starting velocity {vX, vY, vZ} in Cartesian space (m/s).
	*/
	std::array<double, 3> get_birth_vXYZ(const Background::Background& bkg, 
		const Options::Options& opts, const double t, const double x, 
		const double y, const double z, const int tidx, const int xidx,
		const int yidx, const int zidx, pcg32& rng)
	{
		// If we're running a test we will want to start at predetermined
		// velocities to make sure things work correctly. Ideally there
		// wouldn't be an if statement to favor vectorization, but seeing as
		// this is only happens once per particle it's not a huge deal.
		if (opts.bkg_source_int() == 0)
		{
			// Give particles an initial Y velocity to kick off gyration
			if (opts.test_opt_int() == 0 || opts.test_opt_int() == 1
				|| opts.test_opt_int() == 2 || opts.test_opt_int() == 3)
			{
				return {0.0, 10000.0, 0.0};
			}

			// Curvature drift test requires an initial velocity parallel to
			// the field line, we set that to be 30,000 m/s. In this geometry
			// x = R, y = Z and z = phi
			else if (opts.test_opt_int() == 4)
			{
				constexpr double v_par = 30000;
				double v_X {- v_par * std::sin(z)};
				double v_Y {v_par * std::cos(z)};
				return {v_X, v_Y, 0.0};
			}

			// Friction force test case we start at rest so it can accelerate 
			// up to the background velocity
			else if (opts.test_opt_int() == 5)
			{
				return {0.0, 0.0, 0.0};
			}
		}

		// Start at input temperature value
		double start_temp {};
		if (opts.imp_temp_start_opt_int() == 0)
		{
			start_temp = opts.imp_temp_start_val();
		}

		// Start at main ion temperature
		else if (opts.imp_temp_start_opt_int() == 1)
		{
			start_temp = bkg.interp_ti(t, x, y, z); 
			start_temp = std::max(start_temp, opts.min_ti());
		}

		// Sample from a Maxwellian with sigma = sqrt(kT/m) and mean = 0 for an
		// isotropic velocity distribution. T [eV], m [kg]
		const double sigma {std::sqrt(start_temp * Constants::ev_to_j 
			/ (opts.imp_mass_amu() * Constants::amu_to_kg))};  // m/s

		double vX {rng.normal(0.0, sigma)};
		double vY {rng.normal(0.0, sigma)};
		double vZ {rng.normal(0.0, sigma)};

		return {vX, vY, vZ};
	}


	/**
	* Creates a new particle at its starting time, location and velocity.
	*
	* The starting t, x, y and z come from get_birth_val using the
	* imp_*start_opt options. The time index and grid cell indices are then
	* found, the same as in impurity_transport.cpp, and the Cartesian
	* velocity comes from get_birth_vXYZ. The computational velocity is set
	* to zero, since it gets set on the first call to Boris::update_velocity.
	* The weight is 1 / imp_num and the charge is imp_init_charge.
	*
	* @param bkg Background plasma.
	* @param opts Simulation options.
	* @param rng Random number generator for this thread.
	* @return The new particle.
	*/
	ParticleInit make_new_particle(const Background::Background& bkg,
		const Options::Options& opts, pcg32& rng) 
	{
		ParticleInit p;

		// Variables for each grid to avoid repeatedly calling in the loop
		const auto& times = bkg.get_times();
		const auto& grid_x = bkg.get_grid_x();
		const auto& grid_y = bkg.get_grid_y();
		const auto& grid_z = bkg.get_grid_z();

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

		// Get indices, same process as in impurity_transport.cpp
		p.tidx = Indices::get_nearest_index_cpu(times, p.t);
		p.xidx = Indices::get_nearest_cell_index_cpu(grid_x, p.x);
		p.yidx = Indices::get_nearest_cell_index_cpu(grid_y, p.y);
		p.zidx = Indices::get_nearest_cell_index_cpu(grid_z, p.z);

		// Get velocity components in Cartesian space
		auto [vX, vY, vZ] = get_birth_vXYZ(bkg, opts, p.t, p.x, p.y, p.z, 
			p.tidx, p.xidx, p.yidx, p.zidx, rng);
		p.vX = vX;
		p.vY = vY;
		p.vZ = vZ;

		// These get set on the first call to Boris::update_velocity 
		p.vx = 0.0;
		p.vy = 0.0;
		p.vz = 0.0;
		
		// Particle weights are in units of particles/s, meant to represent a
		// given source of particles at this location. We need to normalize
		// by number of particle launched as well.
		p.weight = 1.0 / opts.imp_num();
		p.q = opts.imp_init_charge(); 

		return p;
	}

	/**
	* Replaces dead particles with new ones, as long as there are particles
	* remaining to be followed.
	*
	* The dead slots are counted first. The number refilled is the smaller
	* of the dead count and rem_parts, and rem_parts is reduced by that
	* number. Each refilled slot gets a new particle from make_new_particle,
	* using the random number generator of the thread doing the work, with
	* its particle time reset to 0 and its state set to alive.
	*
	* @param slots Particle slots to refill.
	* @param rem_parts Number of particles still to be followed. Reduced by
	*        the number of particles created, and never goes below 0.
	* @param alive_slots Set to the number of alive slots after refilling.
	* @param bkg Background plasma.
	* @param opts Simulation options.
	* @param rngs Random number generators, one per OpenMP thread.
	*/
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
					slots.set_pt(i, 0.0);  // particle time
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
					slots.set_weight(i, p.weight);
					slots.set_state(i, 0);
				}
			}  // i loop
		}  // pragma omp parallel

		if (rem_parts < 0) rem_parts = 0;

		// Number of alive slots is the number of zeros
		alive_slots = std::count(slots.state().begin(), slots.state().end(), 0);
	}

}  // namespace Slots
