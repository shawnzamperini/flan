/**
* @file  impurity_stats.cpp
*
* @brief Class related to handling impurity transport statistics.
*
* During the simulation, each alive particle adds its count, weight * dt,
* charge and velocities (each multiplied by weight * dt) to the cell and time
* frame it is in. At the end, the statistics from every thread, GPU and MPI
* rank are summed and turned into the impurity density and weighted average
* velocity, charge and Nanbu collisionality in each cell.
*/

#include <cmath>
#include <iomanip>
#include <numeric>
#include <omp.h>
#include <vector>

#include "background.h"
#include "constants.h"
#include "flan_types.h"
#include "impurity_stats.h"
#include "impurity_stats_device.h"
#include "options.h"
#include "slots.h"
#include "vectors.h"

#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

// To-do: Change this namespace to Statistics
namespace Impurity
{
	/**
	* @class Statistics
	*
	* Class to hold all the vectors related to tracking the statistics of
	* an impurity transport simulation. Each array is indexed by
	* (time frame, x cell, y cell, z cell).
	*/

	/**
	* Creates an empty Statistics object with no allocated arrays.
	*/
	Statistics::Statistics() {};

	/**
	* Creates a Statistics object with zero-initialized arrays of size
	* dim1 x dim2 x dim3 x dim4.
	*
	* All the arrays are allocated except the density, which is only
	* allocated in calc_density so it doesn't take up memory for the whole
	* simulation.
	*
	* @param dim1 Number of time frames.
	* @param dim2 Number of grid cells in x.
	* @param dim3 Number of grid cells in y.
	* @param dim4 Number of grid cells in z.
	*/
	Statistics::Statistics(const int dim1, const int dim2, const int dim3, 
		const int dim4)
		: m_dim1 {dim1}
		, m_dim2 {dim2}
		, m_dim3 {dim3}
		, m_dim4 {dim4}
		, m_counts (Vectors::Vector4D<int> {dim1, dim2, dim3, dim4})
		, m_weights (Vectors::Vector4D<BkgFPType> {dim1, dim2, dim3, dim4})
		//, m_density (Vectors::Vector4D<BkgFPType> {dim1, dim2, dim3, dim4})
		, m_charge (Vectors::Vector4D<BkgFPType> {dim1, dim2, dim3, dim4})
	{
		// These arrays can take up a lot of space. Could consider an option
		// to only allocate and track these if requested.

		// Cartesian velocity components
		m_vX.move_into_data(Vectors::Vector4D<BkgFPType> {dim1, dim2, dim3, 
			dim4});
		m_vY.move_into_data(Vectors::Vector4D<BkgFPType> {dim1, dim2, dim3, 
			dim4});
		m_vZ.move_into_data(Vectors::Vector4D<BkgFPType> {dim1, dim2, dim3, 
			dim4});

		// Curvilinear velocity components
		m_vx.move_into_data(Vectors::Vector4D<BkgFPType> {dim1, dim2, dim3, 
			dim4});
		m_vy.move_into_data(Vectors::Vector4D<BkgFPType> {dim1, dim2, dim3, 
			dim4});
		m_vz.move_into_data(Vectors::Vector4D<BkgFPType> {dim1, dim2, dim3, 
			dim4});

		// Collisionality strength from Nanbu model
		m_s.move_into_data(Vectors::Vector4D<BkgFPType> {dim1, dim2, dim3, 
			dim4});
	}

	/**
	* Returns the size of the first dimension of the arrays, which is the
	* number of time frames.
	*
	* @return Size of dimension 1.
	*/
	const int Statistics::get_dim1() const {return m_dim1;}

	/**
	* Returns the size of the second dimension of the arrays, which is the
	* number of grid cells in x.
	*
	* @return Size of dimension 2.
	*/
	const int Statistics::get_dim2() const {return m_dim2;}

	/**
	* Returns the size of the third dimension of the arrays, which is the
	* number of grid cells in y.
	*
	* @return Size of dimension 3.
	*/
	const int Statistics::get_dim3() const {return m_dim3;}

	/**
	* Returns the size of the fourth dimension of the arrays, which is the
	* number of grid cells in z.
	*
	* @return Size of dimension 4.
	*/
	const int Statistics::get_dim4() const {return m_dim4;}

	/**
	* Returns the time (s) along the recorded particle track.
	*
	* Only a single particle track is recorded.
	*
	* @return Reference to the track values, one per recorded step.
	*/
	std::vector<float>& Statistics::get_track_t() {return m_track_t;}

	/**
	* Returns the computational x coordinate along the recorded particle track.
	*
	* Only a single particle track is recorded.
	*
	* @return Reference to the track values, one per recorded step.
	*/
	std::vector<float>& Statistics::get_track_x() {return m_track_x;}

	/**
	* Returns the computational y coordinate along the recorded particle track.
	*
	* Only a single particle track is recorded.
	*
	* @return Reference to the track values, one per recorded step.
	*/
	std::vector<float>& Statistics::get_track_y() {return m_track_y;}

	/**
	* Returns the computational z coordinate along the recorded particle track.
	*
	* Only a single particle track is recorded.
	*
	* @return Reference to the track values, one per recorded step.
	*/
	std::vector<float>& Statistics::get_track_z() {return m_track_z;}

	/**
	* Returns the velocity along computational x along the recorded particle
	* track.
	*
	* Only a single particle track is recorded.
	*
	* @return Reference to the track values, one per recorded step.
	*/
	std::vector<float>& Statistics::get_track_vx() {return m_track_vx;}

	/**
	* Returns the velocity along computational y along the recorded particle
	* track.
	*
	* Only a single particle track is recorded.
	*
	* @return Reference to the track values, one per recorded step.
	*/
	std::vector<float>& Statistics::get_track_vy() {return m_track_vy;}

	/**
	* Returns the velocity along computational z along the recorded particle
	* track.
	*
	* Only a single particle track is recorded.
	*
	* @return Reference to the track values, one per recorded step.
	*/
	std::vector<float>& Statistics::get_track_vz() {return m_track_vz;}

	/**
	* Returns the number of times a particle was recorded in each cell and
	* time frame. Incremented with add_counts.
	*
	* @return Reference to the counts array.
	*/
	Vectors::Vector4D<int>& Statistics::get_counts() {return m_counts;};

	/**
	* Returns the summed weight * dt in each cell and time frame, which is
	* the occupancy in number of particles. Incremented with add_weights.
	*
	* @return Reference to the weights array.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_weights() {return m_weights;}

	/**
	* Returns the impurity density in each cell and time frame (m^-3).
	*
	* Filled by calc_density, and unallocated until then.
	*
	* @return Reference to the density array.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_density() {return m_density;}

	/**
	* Returns the Cartesian X velocity in each cell and time frame.
	*
	* During the simulation this is the running sum of velocity * weight *
	* dt. calc_vels turns it into the weighted average velocity.
	*
	* @return Reference to the vX array.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_vX() {return m_vX;}

	/**
	* Returns the Cartesian Y velocity in each cell and time frame.
	*
	* During the simulation this is the running sum of velocity * weight *
	* dt. calc_vels turns it into the weighted average velocity.
	*
	* @return Reference to the vY array.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_vY() {return m_vY;}

	/**
	* Returns the Cartesian Z velocity in each cell and time frame.
	*
	* During the simulation this is the running sum of velocity * weight *
	* dt. calc_vels turns it into the weighted average velocity.
	*
	* @return Reference to the vZ array.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_vZ() {return m_vZ;}

	/**
	* Returns the velocity along computational x in each cell and time frame.
	*
	* During the simulation this is the running sum of velocity * weight *
	* dt. calc_vels turns it into the weighted average velocity.
	*
	* @return Reference to the vx array.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_vx() {return m_vx;}

	/**
	* Returns the velocity along computational y in each cell and time frame.
	*
	* During the simulation this is the running sum of velocity * weight *
	* dt. calc_vels turns it into the weighted average velocity.
	*
	* @return Reference to the vy array.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_vy() {return m_vy;}

	/**
	* Returns the velocity along computational z in each cell and time frame.
	*
	* During the simulation this is the running sum of velocity * weight *
	* dt. calc_vels turns it into the weighted average velocity.
	*
	* @return Reference to the vz array.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_vz() {return m_vz;}

	/**
	* Returns the charge state in each cell and time frame.
	*
	* During the simulation this is the running sum of charge * weight * dt.
	* calc_charge turns it into the weighted average charge.
	*
	* @return Reference to the charge array.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_charge() {return m_charge;}

	/**
	* Returns the Nanbu collisionality s in each cell and time frame.
	*
	* During the simulation this is the running sum of s * weight * dt.
	* calc_s turns it into the weighted average s.
	*
	* @return Reference to the s array.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_s() {return m_s;}

	/**
	* Adds two Statistics objects together.
	*
	* The counts, weights, charge, s and velocity arrays are summed. Right
	* now only a single particle track is recorded, so the track is copied
	* from other instead of combined.
	*
	* @param other Statistics object to add to this one. Must have the same
	*        dimensions.
	* @return New Statistics object holding the sums.
	*/
	Statistics Statistics::operator+ (const Statistics& other) const
	{
		Statistics ret_stats {m_dim1, m_dim2, m_dim3, m_dim4};
		ret_stats.m_counts = m_counts + other.m_counts;
		ret_stats.m_weights = m_weights + other.m_weights;
		ret_stats.m_charge = m_charge + other.m_charge;
		ret_stats.m_s = m_s + other.m_s;
		ret_stats.m_vX = m_vX + other.m_vX;
		ret_stats.m_vY = m_vY + other.m_vY;
		ret_stats.m_vZ = m_vZ + other.m_vZ;
		ret_stats.m_vx = m_vx + other.m_vx;
		ret_stats.m_vy = m_vy + other.m_vy;
		ret_stats.m_vz = m_vz + other.m_vz;

		// Right now we only track a single particle track, so we don't actually
		// do anything here. This will mean the arrays are empty if there is
		// more than one thread.
		ret_stats.m_track_t = other.m_track_t;
		ret_stats.m_track_x = other.m_track_x;
		ret_stats.m_track_y = other.m_track_y;
		ret_stats.m_track_z = other.m_track_z;
		ret_stats.m_track_vx = other.m_track_vx;
		ret_stats.m_track_vy = other.m_track_vy;
		ret_stats.m_track_vz = other.m_track_vz;

		return ret_stats;
	}

	/**
	* Adds to the counts in a cell and time frame.
	*
	* @param tidx Index of the time frame.
	* @param xidx x index of the grid cell.
	* @param yidx y index of the grid cell.
	* @param zidx z index of the grid cell.
	* @param value Number to add, usually just 1.
	*/
	void Statistics::add_counts(const int tidx, const int xidx, 
		const int yidx, 
		const int zidx, const int value)
	{
		m_counts(tidx, xidx, yidx, zidx) += value;
	}
	
	/**
	* Adds to the weights in a cell and time frame.
	*
	* @param tidx Index of the time frame.
	* @param xidx x index of the grid cell.
	* @param yidx y index of the grid cell.
	* @param zidx z index of the grid cell.
	* @param value Weight to add, usually the particle weight * dt.
	*/
	void Statistics::add_weights(const int tidx, const int xidx, 
		const int yidx, const int zidx, const BkgFPType value)
	{
		m_weights(tidx, xidx, yidx, zidx) += value;
	}

	/**
	* Adds velocities to the running sums in a cell and time frame.
	*
	* @param tidx Index of the time frame.
	* @param xidx x index of the grid cell.
	* @param yidx y index of the grid cell.
	* @param zidx z index of the grid cell.
	* @param vX Cartesian X velocity to add (m/s).
	* @param vY Cartesian Y velocity to add (m/s).
	* @param vZ Cartesian Z velocity to add (m/s).
	* @param vx Velocity along computational x to add.
	* @param vy Velocity along computational y to add.
	* @param vz Velocity along computational z to add.
	* @param bkg Background plasma. Currently unused.
	*/
	void Statistics::add_vels(const int tidx, const int xidx, const int yidx,
		const int zidx, const BkgFPType vX, const BkgFPType vY, 
		const BkgFPType vZ, const BkgFPType vx, const BkgFPType vy, 
		const BkgFPType vz,const Background::Background& bkg)
	{
		// Add velocities to the running total at each cell location
		m_vX(tidx, xidx, yidx, zidx) += vX;
		m_vY(tidx, xidx, yidx, zidx) += vY;
		m_vZ(tidx, xidx, yidx, zidx) += vZ;
		m_vx(tidx, xidx, yidx, zidx) += vx;
		m_vy(tidx, xidx, yidx, zidx) += vy;
		m_vz(tidx, xidx, yidx, zidx) += vz;

		// Note that one could add velocity distribution counting for a future
		// update if the memory load isn't too demanding, but it may be.
	}

	/**
	* Adds to the charge in a cell and time frame.
	*
	* @param tidx Index of the time frame.
	* @param xidx x index of the grid cell.
	* @param yidx y index of the grid cell.
	* @param zidx z index of the grid cell.
	* @param value Charge to add, usually the charge * weight * dt.
	*/
	void Statistics::add_charge(const int tidx, const int xidx, 
		const int yidx, const int zidx, const BkgFPType value)
	{
		m_charge(tidx, xidx, yidx, zidx) += value;
	}

	/**
	* Adds to the Nanbu collisionality s in a cell and time frame.
	*
	* @param tidx Index of the time frame.
	* @param xidx x index of the grid cell.
	* @param yidx y index of the grid cell.
	* @param zidx z index of the grid cell.
	* @param value s to add, usually s * weight * dt.
	*/
	void Statistics::add_s(const int tidx, const int xidx, 
		const int yidx, const int zidx, const BkgFPType value)
	{
		m_s(tidx, xidx, yidx, zidx) += value;
	}

	/**
	* Packs all the arrays into a single buffer that can be sent across MPI
	* processes.
	*
	* The arrays are appended in the order counts, weights, charge, s, vX,
	* vY, vZ, vx, vy, vz. The counts are converted to double. The track and
	* density are not included.
	*
	* @param buf Buffer to fill. Any previous contents are cleared.
	*/
	void Statistics::pack(std::vector<double>& buf) const
	{
		// Clear buffer and reserve enough space for all the Vector4Ds
		buf.clear();
		buf.reserve(m_dim1 * m_dim2 * m_dim3 * m_dim4 * 9); // 9 fields

		// Define a lambda function that breaks down a Vector4D and appends
		// it to our buffer.
		auto append = [&](auto& v4d) {
			for (int i = 0; i < m_dim1; i++)
			for (int j = 0; j < m_dim2; j++)
			for (int k = 0; k < m_dim3; k++)
			for (int l = 0; l < m_dim4; l++)
				buf.push_back(v4d(i,j,k,l));
		};

		// Append each Vector4D to buffer
		append(m_counts);   // convert int → double automatically
		append(m_weights);
		append(m_charge);
		append(m_s);
		append(m_vX);
		append(m_vY);
		append(m_vZ);
		append(m_vx);
		append(m_vy);
		append(m_vz);
	}

	/**
	* Unpacks the arrays from a buffer filled by pack, overwriting the
	* current values.
	*
	* The arrays must already be allocated with the same dimensions as the
	* object that was packed.
	*
	* @param buf Buffer filled by pack.
	*/
	void Statistics::unpack(const std::vector<double>& buf)
	{
		// Index in the 1D buffer
		size_t idx = 0;

		// Define a lambda function to take element in the buffer and put them
		// into a Vector4D.
		auto extract = [&](auto& v4d) {
			for (int i = 0; i < m_dim1; i++)
			for (int j = 0; j < m_dim2; j++)
			for (int k = 0; k < m_dim3; k++)
			for (int l = 0; l < m_dim4; l++)
				v4d(i,j,k,l) = buf[idx++];
		};

		// Extract Vector4D's from buffer. Order must match that in pack!!!
		extract(m_counts);
		extract(m_weights);
		extract(m_charge);
		extract(m_s);
		extract(m_vX);
		extract(m_vY);
		extract(m_vZ);
		extract(m_vx);
		extract(m_vy);
		extract(m_vz);
	}


	/**
	* Calculates the impurity density in each cell and time frame from the
	* weights.
	*
	* The weights were accumulated as weight * dt, so they are already in
	* number of particles, and the particle weights were normalized by the
	* number of particles when each particle was created. The density is
	* then just the weight divided by the cell volume in physical space,
	* J * dx * dy * dz. The density array is allocated here.
	*
	* @param bkg Background plasma, used for the grid and Jacobian.
	* @param tot_imp_num Total number of impurities followed. Currently
	*        unused.
	* @param imp_source_scale_fact Number of particles/s a single particle
	*        represents. Currently unused.
	*/
	void Statistics::calc_density(const Background::Background& bkg, 
		const int tot_imp_num, const double imp_source_scale_fact)
	{
		// Allocate empty density Vector4D. Important we do this here,
		// since if we do it when ImpurityStats gets constructed that means
		// an unused allocated Vector4D is held for the enitre simulation.
		m_density = Vectors::Vector4D<BkgFPType>(m_dim1, m_dim2, m_dim3, 
			m_dim4);

		// Need to loop through the entire Vector4D. Unconventional 
		// indentation here just to keep it clean.
		for (int i {}; i < m_dim1; ++i)
		{
		for (int j {}; j < m_dim2; ++j)
		{
		for (int k {}; k < m_dim3; ++k)
		{
		for (int l {}; l < m_dim4; ++l)
		{
			// Cell width in each direction 
			double dx {bkg.get_grid_x()[j+1] - bkg.get_grid_x()[j]};
			double dy {bkg.get_grid_y()[k+1] - bkg.get_grid_y()[k]};
			double dz {bkg.get_grid_z()[l+1] - bkg.get_grid_z()[l]};

			// Volume in physical space is jacob * dx * dy * dz
			double cell_vol {bkg.get_J()(j,k,l) * dx * dy * dz};

			// Normalize each weight value by the volume and total number of
			// particles launched. This goes from units of (s) to (s/m3). 
			// m_weights is in units of particles/s, so to get total number of
			// particles (weight -> N) we multiply by dt. 
			//m_density(i,j,k,l) = m_weights(i,j,k,l) / cell_vol 
			//	/ tot_imp_num * imp_source_scale_fact;

			// Weights ([particles] since it was accumulated as p_w * dt), so 
			// density is just weight / volume. The weight was already 
			// normalized by the number of particles when the particle was
			// created.
			m_density(i,j,k,l) = m_weights(i,j,k,l) / cell_vol;
		}
		}
		}	
		}

	}
	
	/**
	* Converts the summed velocities in each cell and time frame
	* into weighted average velocities.
	*
	* Each sum is divided by the weight in the cell. Cells with no counts
	* are set to 0. 
	*/
	void Statistics::calc_vels()
	{
		// Need to loop through the entire Vector4D. Unconventional 
		// indentation here just to keep it clean.
		for (int i {}; i < m_dim1; ++i)
		{
		for (int j {}; j < m_dim2; ++j)
		{
		for (int k {}; k < m_dim3; ++k)
		{
		for (int l {}; l < m_dim4; ++l)
		{
			// Total weight (actually, occupancy since this is w * dt)
			double tot_weight {m_weights(i,j,k,l)};

			// Convert each velocity sum in each cell to an average velocity
			// by dividing it by the weight
			double counts {static_cast<double>(m_counts(i,j,k,l))};
			if (counts > 0)
			{
				m_vX(i,j,k,l) /= tot_weight;
				m_vY(i,j,k,l) /= tot_weight;
				m_vZ(i,j,k,l) /= tot_weight;
				m_vx(i,j,k,l) /= tot_weight;
				m_vy(i,j,k,l) /= tot_weight;
				m_vz(i,j,k,l) /= tot_weight;
				
			}
			else
			{
				m_vX(i,j,k,l) = 0.0;
				m_vY(i,j,k,l) = 0.0;
				m_vZ(i,j,k,l) = 0.0;
				m_vx(i,j,k,l) = 0.0;
				m_vy(i,j,k,l) = 0.0;
				m_vz(i,j,k,l) = 0.0;
			}
		}
		}
		}	
		}
	}

	/**
	* Converts the summed charge in each cell and time frame into the
	* weighted average charge.
	*
	* Each sum is divided by the weight in the cell. Cells with no counts are
	* set to 0.
	*/
	void Statistics::calc_charge()
	{

		// Average charge is just the running sum divided by the number of
		// counts in the cell.
		for (int i {}; i < m_dim1; ++i)
		{
		for (int j {}; j < m_dim2; ++j)
		{
		for (int k {}; k < m_dim3; ++k)
		{
		for (int l {}; l < m_dim4; ++l)
		{

			// Total weight (actually, occupancy since this is w * dt)
			double tot_weight {m_weights(i,j,k,l)};

			int counts {m_counts(i,j,k,l)};
			if (counts > 0)
			{
				m_charge(i,j,k,l) /= tot_weight;
			}
			else
			{
				m_charge(i,j,k,l) = 0.0;
			}
		}
		}
		}	
		}
	}

	/**
	* Converts the summed Nanbu collisionality s in each cell and time frame
	* into the weighted average s.
	*
	* Each sum is divided by the weight in the cell. Cells with no counts are
	* set to 0.
	*/
	void Statistics::calc_s()
	{
		// Average s is a weighted average. m_s has accumulated s*w
		for (int i {}; i < m_dim1; ++i)
		{
		for (int j {}; j < m_dim2; ++j)
		{
		for (int k {}; k < m_dim3; ++k)
		{
		for (int l {}; l < m_dim4; ++l)
		{

			// Total weight (actually, occupancy since this is w * dt)
			double tot_weight {m_weights(i,j,k,l)};
			int counts {m_counts(i,j,k,l)};
			if (counts > 0)
			{
				m_s(i,j,k,l) /= tot_weight;
			}
			else
			{
				m_s(i,j,k,l) = 0.0;
			}
		}
		}
		}	
		}
	}


	/**
	* Allocates GPU memory for the statistics arrays and copies them to the
	* device.
	*
	* The counts, weights, velocity, charge and s arrays are copied. The
	* track and density are not. If GPU support was not compiled in, an
	* empty StatisticsDevice is returned.
	*
	* @param device_id ID of the GPU to allocate the memory on.
	* @return StatisticsDevice struct holding the device pointers. The memory
	*         must be released with free_stats.
	*/
	ImpurityStats::StatisticsDevice Statistics::to_device(int device_id)
	{

		// Create empty device-side structure
		ImpurityStats::StatisticsDevice stats_d {};

#ifdef USE_CUDA

		stats_d.device_id = device_id;
		stats_d.Nt = m_dim1;
		stats_d.Nx = m_dim2;
		stats_d.Ny = m_dim3;
		stats_d.Nz = m_dim4;
		int N {m_dim1 * m_dim2 * m_dim3 * m_dim4};

		// GPU allocation
		cudaSetDevice(device_id);

		cudaMalloc(&stats_d.counts, N * sizeof(int));
		cudaMalloc(&stats_d.weights, N * sizeof(double));
		cudaMalloc(&stats_d.vX, N * sizeof(double));
		cudaMalloc(&stats_d.vY, N * sizeof(double));
		cudaMalloc(&stats_d.vZ, N * sizeof(double));
		cudaMalloc(&stats_d.vx, N * sizeof(double));
		cudaMalloc(&stats_d.vy, N * sizeof(double));
		cudaMalloc(&stats_d.vz, N * sizeof(double));
		cudaMalloc(&stats_d.q, N * sizeof(double));
		cudaMalloc(&stats_d.s, N * sizeof(double));

		// Copy to device
		cudaMemcpy(stats_d.counts, m_counts.get_data().data(), N * sizeof(int), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(stats_d.weights, m_weights.get_data().data(), N * sizeof(double), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(stats_d.vX, m_vX.get_data().data(), N * sizeof(double), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(stats_d.vY, m_vY.get_data().data(), N * sizeof(double), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(stats_d.vZ, m_vZ.get_data().data(), N * sizeof(double), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(stats_d.vx, m_vx.get_data().data(), N * sizeof(double), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(stats_d.vy, m_vy.get_data().data(), N * sizeof(double), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(stats_d.vz, m_vz.get_data().data(), N * sizeof(double), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(stats_d.q, m_charge.get_data().data(), N * sizeof(double), 
			cudaMemcpyHostToDevice);
		cudaMemcpy(stats_d.s, m_s.get_data().data(), N * sizeof(double), 
			cudaMemcpyHostToDevice);

#endif

		return stats_d;

	}

	/**
	* Copies the statistics from a GPU and adds them to this object.
	*
	* Waits for the device to finish, copies each array back to the host and
	* adds it to the matching array here, printing progress as it goes. Does
	* nothing if GPU support was not compiled in.
	*
	* @param stats_d StatisticsDevice struct returned by to_device.
	* @param device_id ID of the GPU holding the statistics.
	*/
	void Statistics::add_stats_device(ImpurityStats::StatisticsDevice& stats_d,
		int device_id)
	{

#ifdef USE_CUDA

		// Allocate vectors to hold GPU data, these are the middleman
		int N {m_dim1 * m_dim2 * m_dim3 * m_dim4};
		std::vector<int> counts_h (N);
		std::vector<double> weights_h (N);
		std::vector<double> vX_h (N);
		std::vector<double> vY_h (N);
		std::vector<double> vZ_h (N);
		std::vector<double> vx_h (N);
		std::vector<double> vy_h (N);
		std::vector<double> vz_h (N);
		std::vector<double> s_h (N);
		std::vector<double> charge_h (N);

		// Previous errors could cause this to hang, unsure how to handle
		// this yet
		cudaError_t err = cudaDeviceSynchronize();
		printf("sync: %s\n", cudaGetErrorString(err));

		// Copy over from GPU to host
		std::cout << "Copying data from GPU #" << device_id << " (N = " 
			<< N << ")...\n";
		cudaSetDevice(device_id);
		std::cout << "  - counts\n";
		cudaMemcpy(counts_h.data(), stats_d.counts, N * sizeof(int), 
			cudaMemcpyDeviceToHost);
		std::cout << "  - weights\n";
		cudaMemcpy(weights_h.data(), stats_d.weights, N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		std::cout << "  - vX\n";
		cudaMemcpy(vX_h.data(), stats_d.vX, N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		std::cout << "  - vY\n";
		cudaMemcpy(vY_h.data(), stats_d.vY, N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		std::cout << "  - vZ\n";
		cudaMemcpy(vZ_h.data(), stats_d.vZ, N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		std::cout << "  - vx\n";
		cudaMemcpy(vx_h.data(), stats_d.vx, N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		std::cout << "  - vy\n";
		cudaMemcpy(vy_h.data(), stats_d.vy, N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		std::cout << "  - vz\n";
		cudaMemcpy(vz_h.data(), stats_d.vz, N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		std::cout << "  - Nanbu s\n";
		cudaMemcpy(s_h.data(), stats_d.s, N * sizeof(double), 
			cudaMemcpyDeviceToHost);
		std::cout << "  - charge\n";
		cudaMemcpy(charge_h.data(), stats_d.q, N * sizeof(double), 
			cudaMemcpyDeviceToHost);

		// Loop through each dimension
		std::cout << "Adding to host object...\n";
		#pragma omp parallel for
		for (int i=0; i < m_dim1; ++i)
		{
		for (int j {}; j < m_dim2; ++j)
		{
		for (int k {}; k < m_dim3; ++k)
		{
		for (int l {}; l < m_dim4; ++l)
		{
			// 4D to 1D index
			int idx {m_counts.calc_index(i,j,k,l)};

			// Add statistics into this object
			m_counts(i,j,k,l) += counts_h[idx];
			m_weights(i,j,k,l) += weights_h[idx];
			m_vX(i,j,k,l) += vX_h[idx];
			m_vY(i,j,k,l) += vY_h[idx];
			m_vZ(i,j,k,l) += vZ_h[idx];
			m_vx(i,j,k,l) += vx_h[idx];
			m_vy(i,j,k,l) += vy_h[idx];
			m_vz(i,j,k,l) += vz_h[idx];
			m_s(i,j,k,l)  += s_h[idx];
			m_charge(i,j,k,l) += charge_h[idx];
		}
		}
		}
		}

#endif

	}

	/**
	* Sums a Statistics object across all MPI ranks.
	*
	* The arrays are packed into a buffer, summed onto rank 0 with
	* MPI_Reduce and unpacked into a new object. With a single process the
	* object is returned as-is.
	*
	* @param local_stats Statistics recorded on this rank.
	* @return On rank 0, the summed statistics. On other ranks, a
	*         zero-initialized object with the same dimensions.
	*/
	Statistics reduce_stats(const Statistics& local_stats)
	{

		// MPI rank
		int rank {};
		int nprocs {};
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

		// No reduction needed with a single process, return as-is
		if (nprocs == 1) return local_stats;

		// Send and recieve buffers
		std::vector<double> sendbuf;
		std::vector<double> recvbuf;

		// Pack up stats local to this rank and resize the recieve buffer to
		// match
		local_stats.pack(sendbuf);
		recvbuf.resize(sendbuf.size());

		// Call reduce to add all the Vector4Ds in Statistics together, putting
		// the result in recvbuf. 
		MPI_Reduce(sendbuf.data(), recvbuf.data(), sendbuf.size(), 
			MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		// Create object to store the summed results and unpack into it
		// Only rank 0 returns a meaningful object
		if (rank == 0) {
			Statistics global_stats {local_stats.get_dim1(), 
				local_stats.get_dim2(), local_stats.get_dim3(), 
				local_stats.get_dim4()};
			global_stats.unpack(recvbuf);
			return global_stats;
		}

		// Other ranks return a dummy object
		return Statistics(local_stats.get_dim1(),
						  local_stats.get_dim2(),
						  local_stats.get_dim3(),
						  local_stats.get_dim4());
	}

	/**
	* Frees the GPU memory held by a StatisticsDevice struct.
	*
	* Every array allocated in Statistics::to_device is freed and its pointer
	* is set to nullptr to avoid accidental reuse. If GPU support was not
	* compiled in, only the pointers are reset.
	*
	* @param stats_d StatisticsDevice struct whose memory is freed.
	* @param device_id ID of the GPU the memory was allocated on.
	*/
	void free_stats(ImpurityStats::StatisticsDevice& stats_d, int device_id)
	{

#ifdef USE_CUDA

		cudaSetDevice(device_id);

		// Free memory
		cudaFree(stats_d.counts);
		cudaFree(stats_d.weights);
		cudaFree(stats_d.vX);
		cudaFree(stats_d.vY);
		cudaFree(stats_d.vZ);
		cudaFree(stats_d.vx);
		cudaFree(stats_d.vy);
		cudaFree(stats_d.vz);
		cudaFree(stats_d.s);
		cudaFree(stats_d.q);

#endif

		// Set pointers to null to avoid accidental reuse
		stats_d.counts = nullptr;
		stats_d.weights = nullptr;
		stats_d.vX = nullptr;
		stats_d.vY = nullptr;
		stats_d.vZ = nullptr;
		stats_d.vx = nullptr;
		stats_d.vy = nullptr;
		stats_d.vz = nullptr;
		stats_d.s = nullptr;
		stats_d.q = nullptr;
	}

	/*
	// buffered approach
	void record_stats_cpu(Statistics& imp_stats, const Slots::Slots& slots,
		const Options::Options& opts, const double imp_time_step)
	{
		// Number of distinct thread-local arrays we need pointers to.
		constexpr int NUM_FIELDS = 9;

		// Sized to the max threads OpenMP could ever hand us -- allocated once,
		// reused every timestep. Each inner vector<double*> holds one pointer
		// per thread for that field.
		static std::vector<std::vector<double*>> field_ptrs(NUM_FIELDS);
		static std::vector<int*> counts_ptrs;  // counts is int, kept separate
		static std::size_t grid_size {0};
		static int num_threads_used {0};

		#pragma omp parallel
		{
			static thread_local std::vector<int> counts_local;
			static thread_local std::vector<double> weights_local;
			static thread_local std::vector<double> vX_local;
			static thread_local std::vector<double> vY_local;
			static thread_local std::vector<double> vZ_local;
			static thread_local std::vector<double> vx_local;
			static thread_local std::vector<double> vy_local;
			static thread_local std::vector<double> vz_local;
			static thread_local std::vector<double> charge_local;

			const int tid {omp_get_thread_num()};

			// First-touch setup: allocate this thread's local buffers, and
			// register this thread's pointers into the shared lookup arrays.
			// Every thread writes only index [tid] of each pointer array, so
			// no two threads ever touch the same slot -- no lock needed here.
			if (counts_local.empty())
			{
				counts_local.resize(imp_stats.get_counts().get_data().size());
				weights_local.resize(imp_stats.get_weights().get_data().size());
				vX_local.resize(imp_stats.get_vX().get_data().size());
				vY_local.resize(imp_stats.get_vY().get_data().size());
				vZ_local.resize(imp_stats.get_vZ().get_data().size());
				vx_local.resize(imp_stats.get_vx().get_data().size());
				vy_local.resize(imp_stats.get_vy().get_data().size());
				vz_local.resize(imp_stats.get_vz().get_data().size());
				charge_local.resize(imp_stats.get_charge().get_data().size());

				// One thread (any one -- doesn't matter which) sizes the
				// shared pointer/bookkeeping arrays before anyone writes into
				// them, so every thread's write below lands in a valid slot.
				#pragma omp single
				{
					num_threads_used = omp_get_num_threads();
					grid_size = counts_local.size();

					counts_ptrs.resize(num_threads_used);
					for (auto& v : field_ptrs) v.resize(num_threads_used);
				}
				// implicit barrier at the end of "single" -- guarantees the
				// resize above has happened before any thread writes its
				// pointer into these arrays below.

				counts_ptrs[tid]   = counts_local.data();
				field_ptrs[0][tid] = weights_local.data();
				field_ptrs[1][tid] = vX_local.data();
				field_ptrs[2][tid] = vY_local.data();
				field_ptrs[3][tid] = vZ_local.data();
				field_ptrs[4][tid] = vx_local.data();
				field_ptrs[5][tid] = vy_local.data();
				field_ptrs[6][tid] = vz_local.data();
				field_ptrs[7][tid] = charge_local.data();
			}

			// Zero local buffers each timestep
			std::fill(counts_local.begin(), counts_local.end(), 0);
			std::fill(weights_local.begin(), weights_local.end(), 0.0);
			std::fill(vX_local.begin(), vX_local.end(), 0.0);
			std::fill(vY_local.begin(), vY_local.end(), 0.0);
			std::fill(vZ_local.begin(), vZ_local.end(), 0.0);
			std::fill(vx_local.begin(), vx_local.end(), 0.0);
			std::fill(vy_local.begin(), vy_local.end(), 0.0);
			std::fill(vz_local.begin(), vz_local.end(), 0.0);
			std::fill(charge_local.begin(), charge_local.end(), 0.0);

			// Independent loops, so nowait
			#pragma omp for nowait
			for (int i = 0; i < slots.N(); ++i)
			{
				if (slots.state()[i]) continue;

				int idx {imp_stats.get_counts().calc_index(slots.tidx()[i],
					slots.xidx()[i], slots.yidx()[i], slots.zidx()[i])};

				double p_w {slots.weight()[i]};
				counts_local[idx] += 1;
				weights_local[idx] += p_w * imp_time_step;

				vX_local[idx] += slots.vX()[i] * p_w;
				vY_local[idx] += slots.vY()[i] * p_w;
				vZ_local[idx] += slots.vZ()[i] * p_w;
				vx_local[idx] += slots.vx()[i] * p_w;
				vy_local[idx] += slots.vy()[i] * p_w;
				vz_local[idx] += slots.vz()[i] * p_w;
				charge_local[idx] += slots.q()[i] * p_w;
			}

			// Barrier: every thread must finish its "omp for" (including the
			// implicit end-of-for wait would normally happen here, but we used
			// nowait above) before the merge below reads any thread's buffer.
			#pragma omp barrier

			// Merge: parallelized over grid index j instead of serialized over
			// threads. Each thread sums every other thread's contribution for
			// its own slice of j -- so this now gets *cheaper* per thread as
			// thread count grows, rather than more expensive.
			#pragma omp for schedule(static)
			for (std::size_t j = 0; j < grid_size; ++j)
			{
				int    count_sum  {0};
				double weight_sum {0.0}, vX_sum {0.0}, vY_sum {0.0}, vZ_sum {0.0};
				double vx_sum {0.0}, vy_sum {0.0}, vz_sum {0.0}, charge_sum {0.0};

				for (int t = 0; t < num_threads_used; ++t)
				{
					count_sum  += counts_ptrs[t][j];
					weight_sum += field_ptrs[0][t][j];
					vX_sum     += field_ptrs[1][t][j];
					vY_sum     += field_ptrs[2][t][j];
					vZ_sum     += field_ptrs[3][t][j];
					vx_sum     += field_ptrs[4][t][j];
					vy_sum     += field_ptrs[5][t][j];
					vz_sum     += field_ptrs[6][t][j];
					charge_sum += field_ptrs[7][t][j];
				}

				imp_stats.get_counts().get_data()[j]  += count_sum;
				imp_stats.get_weights().get_data()[j] += weight_sum;
				imp_stats.get_vX().get_data()[j]      += vX_sum;
				imp_stats.get_vY().get_data()[j]      += vY_sum;
				imp_stats.get_vZ().get_data()[j]      += vZ_sum;
				imp_stats.get_vx().get_data()[j]      += vx_sum;
				imp_stats.get_vy().get_data()[j]      += vy_sum;
				imp_stats.get_vz().get_data()[j]      += vz_sum;
				imp_stats.get_charge().get_data()[j]  += charge_sum;
			}
			// implicit barrier at end of this "omp for" -- fine, since this is
			// the last thing the parallel region does anyway.

		} // pragma omp parallel

	} // record_stats_cpu
	*/

	/**
	* Records the statistics of every alive particle into its cell and time
	* frame.
	*
	* For each particle, the count is incremented and weight * dt (the
	* occupancy, in number of particles) is added to the weights. The
	* velocities and charge are added multiplied by weight * dt, so they can
	* later be turned into weighted averages.
	*
	* Every thread writes straight into the shared arrays with atomic
	* updates. In my basic testing this approach is much faster than the
	* buffered approach (commented out above), though it's possible the
	* buffered approach wins ground when there are significantly more
	* particles than grid cells. That's because in that case particles may
	* be fighting over accumulating stats in a cell, which can get held up in
	* the atomics here. I am partial to the atomic approach just because it's
	* so much cleaner.
	*
	* @param imp_stats Statistics to record into.
	* @param slots Particle slots.
	* @param opts Simulation options. Currently unused.
	* @param dt Time step (s).
	*/
	void record_stats_cpu(Statistics& imp_stats, const Slots::Slots& slots,
		const Options::Options& opts, const double dt)
	{
		// Grab raw pointers once outside the parallel region -- avoids
		// repeated virtual/accessor call overhead inside the hot loop.
		int*    counts_data  {imp_stats.get_counts().get_data().data()};
		double* weights_data {imp_stats.get_weights().get_data().data()};
		double* vX_data      {imp_stats.get_vX().get_data().data()};
		double* vY_data      {imp_stats.get_vY().get_data().data()};
		double* vZ_data      {imp_stats.get_vZ().get_data().data()};
		double* vx_data      {imp_stats.get_vx().get_data().data()};
		double* vy_data      {imp_stats.get_vy().get_data().data()};
		double* vz_data      {imp_stats.get_vz().get_data().data()};
		double* charge_data  {imp_stats.get_charge().get_data().data()};

		// No thread-local buffers, no zero-fill, no merge pass. Every thread
		// writes straight into the shared arrays; atomics make each individual
		// += race-free. Total work is now O(N particles), not O(grid_size) --
		// grid cells that see no activity this timestep cost nothing.
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < slots.N(); ++i)
		{
			if (slots.state()[i]) continue;

			int idx {imp_stats.get_counts().calc_index(slots.tidx()[i],
				slots.xidx()[i], slots.yidx()[i], slots.zidx()[i])};

			// With weights in [particles/s], this is then number of particles.
			// We can call this the occupancy.
			double p_w_dt {slots.weight()[i] * dt};

			#pragma omp atomic update
			counts_data[idx] += 1;

			#pragma omp atomic update
			weights_data[idx] += p_w_dt;

			#pragma omp atomic update
			vX_data[idx] += slots.vX()[i] * p_w_dt;

			#pragma omp atomic update
			vY_data[idx] += slots.vY()[i] * p_w_dt;

			#pragma omp atomic update
			vZ_data[idx] += slots.vZ()[i] * p_w_dt;

			#pragma omp atomic update
			vx_data[idx] += slots.vx()[i] * p_w_dt;

			#pragma omp atomic update
			vy_data[idx] += slots.vy()[i] * p_w_dt;

			#pragma omp atomic update
			vz_data[idx] += slots.vz()[i] * p_w_dt;

			#pragma omp atomic update
			charge_data[idx] += slots.q()[i] * p_w_dt;

		}

	} // record_stats_cpu

}
