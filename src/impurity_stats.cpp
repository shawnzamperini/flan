/**
* @file impurity_stats.cpp
*
* @brief Class related to handling impurity transport statistics
*/

#include <cmath>
#include <iomanip>
#include <numeric>
#include <vector>

#include "background.h"
#include "constants.h"
#include "flan_types.h"
#include "impurity.h"
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
	* @brief Impurity transport statistics class
	*
	* Class to hold all the vectors related to tracking the statistics of
	* an impurity transport simulation.
	*/

	Statistics::Statistics() {};

	/**
	* @brief Constructor
	*
	* @param dim1 Size of first dimension of containing Vector4D's
	* @param dim2 Size of second dimension of containing Vector4D's
	* @param dim3 Size of third dimension of containing Vector4D's
	* @param dim4 Size of fourth dimension of containing Vector4D's
	*
	* Zero-initializes the member Vector4D's.
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
	* @brief Accessor for dim1
	*/
	const int Statistics::get_dim1() const {return m_dim1;}
	const int Statistics::get_dim2() const {return m_dim2;}
	const int Statistics::get_dim3() const {return m_dim3;}
	const int Statistics::get_dim4() const {return m_dim4;}

	/**
	* @brief Accessor for particle track time
	*/
	std::vector<float>& Statistics::get_track_t() {return m_track_t;}

	/**
	* @brief Accessor for particle track x position
	*/
	std::vector<float>& Statistics::get_track_x() {return m_track_x;}

	/**
	* @brief Accessor for particle track y position
	*/
	std::vector<float>& Statistics::get_track_y() {return m_track_y;}

	/**
	* @brief Accessor for particle track z position
	*/
	std::vector<float>& Statistics::get_track_z() {return m_track_z;}

	/**
	* @brief Accessor for particle track x velocity
	*/
	std::vector<float>& Statistics::get_track_vx() {return m_track_vx;}

	/**
	* @brief Accessor for particle track y velocity
	*/
	std::vector<float>& Statistics::get_track_vy() {return m_track_vy;}

	/**
	* @brief Accessor for particle track z velocity
	*/
	std::vector<float>& Statistics::get_track_vz() {return m_track_vz;}

	/**
	* @brief Accessor for counts data
	* @return Vector4D<int>
	* @sa add_counts()
	*
	* add_counts() is used to increment data in this vector.
	*/
	Vectors::Vector4D<int>& Statistics::get_counts() {return m_counts;};

	/**
	* @brief Accessor for weights data
	* @return Vector4D<BkgFPType>
	* @sa add_weights()
	*
	* add_weights() is used to add data into this vector
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_weights() {return m_weights;}

	/**
	* @brief Accessor for density data
	* @return Vector4D<BkgFPType>
	* @sa calc_density()
	*
	* calc_density() is used to fill this vector. Remains unallocated until
	* then.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_density() {return m_density;}

	/**
	* @brief Accessor for X velocity data
	* @return Vector4D<BkgFPType>
	* @sa calc_vels()
	*
	* calc_vels() is used to turn the data in this into an actual velocity.
	* Before doing that it's just a collection of Monte Carlo type values.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_vX() {return m_vX;}

	/**
	* @brief Accessor for Y velocity data
	* @return Vector4D<BkgFPType>
	* @sa calc_vels()
	*
	* calc_vels() is used to turn the data in this into an actual velocity.
	* Before doing that it's just a collection of Monte Carlo type values.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_vY() {return m_vY;}

	/**
	* @brief Accessor for Z velocity data
	* @return Vector4D<BkgFPType>
	* @sa calc_vels()
	*
	* calc_vels() is used to turn the data in this into an actual velocity.
	* Before doing that it's just a collection of Monte Carlo type values.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_vZ() {return m_vZ;}

	/**
	* @brief Accessor for x velocity data
	* @return Vector4D<BkgFPType>
	* @sa calc_vels()
	*
	* calc_vels() is used to turn the data in this into an actual velocity.
	* Before doing that it's just a collection of Monte Carlo type values.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_vx() {return m_vx;}

	/**
	* @brief Accessor for y velocity data
	* @return Vector4D<BkgFPType>
	* @sa calc_vels()
	*
	* calc_vels() is used to turn the data in this into an actual velocity.
	* Before doing that it's just a collection of Monte Carlo type values.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_vy() {return m_vy;}

	/**
	* @brief Accessor for z velocity data
	* @return Vector4D<BkgFPType>
	* @sa calc_vels()
	*
	* calc_vels() is used to turn the data in this into an actual velocity.
	* Before doing that it's just a collection of Monte Carlo type values.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_vz() {return m_vz;}

	/**
	* @brief Accessor for charge data
	* @return Vector4D<BkgFPType>
	* @sa calc_charge()
	*
	* calc_charge() is used to turn the aggregated data into average
	* charge values at each location.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_charge() {return m_charge;}

	/**
	* @brief Accessor for Nanbu collisionality data
	* @return Vector4D<BkgFPType>
	* @sa calc_s()
	*
	* calc_s() is used to turn the aggregated data into average
	* s values at each location.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_s() {return m_s;}

	/**
	* @brief Overload + to add two Statistics together
	* @param other A Statistics object
	* @return Statistics object with the summed data
	*
	* The values of counts and weights are added together and returned in
	* a new Statistics object.
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
	* @brief Increment counts at the specified indices
	* @param tidx Time index
	* @param xidx x index
	* @param yidx y index
	* @param zidx z index
	* @param value Number to increment by (usually just 1)
	*/
	void Statistics::add_counts(const int tidx, const int xidx, 
		const int yidx, 
		const int zidx, const int value)
	{
		m_counts(tidx, xidx, yidx, zidx) += value;
	}
	
	/**
	* @brief Increment weights at the specified indices
	* @param tidx Time index
	* @param xidx x index
	* @param yidx y index
	* @param zidx z index
	* @param value Weight to add to the cell
	*/
	void Statistics::add_weights(const int tidx, const int xidx, 
		const int yidx, const int zidx, const BkgFPType value)
	{
		m_weights(tidx, xidx, yidx, zidx) += value;
	}

	/**
	*
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
	* @brief Increment charge at the specified indices
	* @param tidx Time index
	* @param xidx x index
	* @param yidx y index
	* @param zidx z index
	* @param value Charge to add to the cell
	*/
	void Statistics::add_charge(const int tidx, const int xidx, 
		const int yidx, const int zidx, const BkgFPType value)
	{
		m_charge(tidx, xidx, yidx, zidx) += value;
	}

	/**
	* @brief Increment Nanbu collisionality strength at the specified indices
	* @param tidx Time index
	* @param xidx x index
	* @param yidx y index
	* @param zidx z index
	* @param value s to add to the cell
	*/
	void Statistics::add_s(const int tidx, const int xidx, 
		const int yidx, const int zidx, const BkgFPType value)
	{
		m_s(tidx, xidx, yidx, zidx) += value;
	}

	// Add impurity position to track vectors. This is mainly for testing
	// that the path of a particle follows a particular route and showing that
	// the physics is correct.
	void Statistics::update_track(const Impurity& imp)
	{
		m_track_t.push_back(imp.get_t());
		m_track_x.push_back(imp.get_x());
		m_track_y.push_back(imp.get_y());
		m_track_z.push_back(imp.get_z());
		m_track_vx.push_back(imp.get_vx());
		m_track_vy.push_back(imp.get_vy());
		m_track_vz.push_back(imp.get_vz());
	}



	/**
	* @brief Pack all the arrays into a single buffer that can be sent across
	*        MPI processes
	*
	* @param buf The buffer that will contain all the Vector4Ds in Statistics
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
	* @brief Unpack all the arrays that were packed via pack
	*
	* @param buf The buffer that contains all the Vector4Ds in Statistics
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
	* @brief Calculate impurity density from data stored in counts and weights
	* @param bkg The Background object used in the simulation
	* @param tot_imp_num The total number of impurity ions followed (primary
	* and secondary)
	* @param imp_source_scale_fact A scaling factor related to the impurity
	* source in units of particles/s.\ Essentially says how many particles a 
	* single particle represents.
	*
	* Calculate the density in each cell using the Monte Carlo counts and
	* weights. The Background is needed to get the cell volumes. As of now,
	* imp_source_scale_fact is user supplied, but it hasn't been fully thought
	* out yet. This only works for uniform Cartesian grids right now, but it
	* can be generalized.
	*/
	void Statistics::calc_density(const Background::Background& bkg, 
		const int tot_imp_num, const double imp_source_scale_fact)
	{
		std::cout << "calc_density: NEED TO THINK ABOUT IF THIS IS RIGHT\n";
		// Allocate empty density Vector4D. Important we do this here,
		// since if we do it when ImpurityStats gets constructed that means
		// a pointless allocated Vector4D gets created for each OpenMP thread
		// during the impurity following code!
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
			m_density(i,j,k,l) = m_weights(i,j,k,l) / cell_vol 
				/ tot_imp_num * imp_source_scale_fact;
		}
		}
		}	
		}

	}
	
	/**
	* @brief Calculate the average velocity in each cell
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
			// Total weight
			double tot_weight {m_weights(i,j,k,l)};

			// Convert each velocity sum in each cell to an average velocity
			// by dividing it by the number of counts.
			double counts {static_cast<double>(m_counts(i,j,k,l))};
			if (counts > 0)
			{
				// This is actually wrong for a weighted average
				//m_vX(i,j,k,l) /= counts;
				//m_vY(i,j,k,l) /= counts;
				//m_vZ(i,j,k,l) /= counts;
				m_vX(i,j,k,l) /= tot_weight;
				m_vY(i,j,k,l) /= tot_weight;
				m_vZ(i,j,k,l) /= tot_weight;
				
			}
			else
			{
				m_vX(i,j,k,l) = 0.0;
				m_vY(i,j,k,l) = 0.0;
				m_vZ(i,j,k,l) = 0.0;
			}
		}
		}
		}	
		}
	}

	/**
	* @brief Calculate the average charge at each cell location
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
			int counts {m_counts(i,j,k,l)};
			if (counts > 0)
			{
				m_charge(i,j,k,l) /= counts;
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
	* @brief Calculate the average s at each cell location
	*/
	void Statistics::calc_s()
	{
		// Average s is just the running sum divided by the number of
		// counts in the cell.
		for (int i {}; i < m_dim1; ++i)
		{
		for (int j {}; j < m_dim2; ++j)
		{
		for (int k {}; k < m_dim3; ++k)
		{
		for (int l {}; l < m_dim4; ++l)
		{
			int counts {m_counts(i,j,k,l)};
			if (counts > 0)
			{
				m_s(i,j,k,l) /= counts;
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
		//cudaMalloc(&stats_d.s, N * sizeof(double));

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
		//cudaMemcpy(stats_d.s, m_s.get_data().data(), N * sizeof(double), 
		//	cudaMemcpyHostToDevice);

#endif

		return stats_d;

	}

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
		std::cout << "  - charge\n";
		cudaMemcpy(charge_h.data(), stats_d.q, N * sizeof(double), 
			cudaMemcpyDeviceToHost);

		// Loop through each dimension
		std::cout << "Adding to host object...\n";
		#pragma omp for
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
			m_charge(i,j,k,l) += charge_h[idx];
		}
		}
		}
		}

#endif

	}

	/**
	* @brief Helper function to reduce a Statistic object across MPI ranks
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
		stats_d.q = nullptr;
	}

	void record_stats_cpu(Statistics& imp_stats, const Slots::Slots& slots,
		const Options::Options& opts, const double imp_time_step)
	{

		// We are accumulating statistics in a simple imp_stats object, so
		// we have to be sensitive to race conditions. The solution is to use
		// thread local arrays and accumulate the stats into those each time
		// step, and then add those into the imp_stats array under a critical
		// block. 
		#pragma omp parallel
		{
			// Allocate local arrays for each thread. static tells us they're
			// global and will stay for the lifetime of the program, so we 
			// don't need to waste time reallocating these arrays every time
			// step. thread_local allocates it only for that thread.
			static thread_local std::vector<int> counts_local;
			static thread_local std::vector<double> weights_local;
			static thread_local std::vector<double> vX_local;
			static thread_local std::vector<double> vY_local;
			static thread_local std::vector<double> vZ_local;
			static thread_local std::vector<double> vx_local;
			static thread_local std::vector<double> vy_local;
			static thread_local std::vector<double> vz_local;
			static thread_local std::vector<double> charge_local;

			// Allocate once per thread. This only executes when the thread
			// is created the first time, so only once.
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
			}

			// Zero them at each timestep. This sounds expensive but the 
			// compiler should make this a pretty fast operation. We could
			// add a timer here if we're curious.
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
				// Skip dead particles (0=alive, 1+=dead)
				if (slots.state()[i]) continue;

				// 4D -> 1D index
				int idx {imp_stats.get_counts().calc_index(slots.tidx()[i], 
					slots.xidx()[i], slots.yidx()[i], slots.zidx()[i])};

				// Increment local arrays
				double p_w {slots.weight()[i]};
				counts_local[idx] += 1;
				weights_local[idx] += p_w * imp_time_step;
				
				// Add property * weight since we are doing a weighted average
				vX_local[idx] += slots.vX()[i] * p_w;
				vY_local[idx] += slots.vY()[i] * p_w;
				vZ_local[idx] += slots.vZ()[i] * p_w;
				vx_local[idx] += slots.vx()[i] * p_w;
				vy_local[idx] += slots.vy()[i] * p_w;
				vz_local[idx] += slots.vz()[i] * p_w;
				charge_local[idx] += slots.q()[i] * p_w;
			}

			// Merge into global stats
			#pragma omp critical
			{
				for (std::size_t j = 0; j < counts_local.size(); ++j)
				{
					imp_stats.get_counts().get_data()[j] += counts_local[j];
					imp_stats.get_weights().get_data()[j] += weights_local[j];
					imp_stats.get_vX().get_data()[j] += vX_local[j];
					imp_stats.get_vY().get_data()[j] += vY_local[j];
					imp_stats.get_vZ().get_data()[j] += vZ_local[j];
					imp_stats.get_vx().get_data()[j] += vx_local[j];
					imp_stats.get_vy().get_data()[j] += vy_local[j];
					imp_stats.get_vz().get_data()[j] += vz_local[j];
				}
			}
		} // pragma omp parallel

	} // record_stats_cpu

}
