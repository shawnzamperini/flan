/**
* @file impurity_stats.cpp
*
* @brief Class related to handling impurity transport statistics
*/

#include <numeric>
#include <vector>

#include "vectors.h"
#include "impurity_stats.h"

namespace Impurity
{
	/** 
	* @class Statistics
	* @brief Impurity transport statistics class
	*
	* Class to hold all the vectors related to tracking the statistics of
	* an impurity transport simulation.
	*/

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
		const int dim4, const bool vel_stats)
		: m_dim1 {dim1}
		, m_dim2 {dim2}
		, m_dim3 {dim3}
		, m_dim4 {dim4}
		, m_counts (Vectors::Vector4D<int> {dim1, dim2, dim3, dim4})
		, m_weights (Vectors::Vector4D<double> {dim1, dim2, dim3, dim4})
		, m_density (Vectors::Vector4D<double> {dim1, dim2, dim3, dim4})
		, m_vel_stats (vel_stats)
	{
		// These aren't always needed, can cut back on memory usage by not
		// including them by default.
		if (m_vel_stats)
		{
			m_vx = (Vectors::Vector4D<double> {dim1, dim2, dim3, dim4});
			m_vy = (Vectors::Vector4D<double> {dim1, dim2, dim3, dim4});
			m_vz = (Vectors::Vector4D<double> {dim1, dim2, dim3, dim4});
		}
	}

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
	* @return Vector4D<double>
	* @sa add_weights()
	*
	* add_weights() is used to add data into this vector
	*/
	Vectors::Vector4D<double>& Statistics::get_weights() {return m_weights;}

	/**
	* @brief Accessor for density data
	* @return Vector4D<double>
	* @sa calc_density()
	*
	* calc_density() is used to fill this vector.
	*/
	Vectors::Vector4D<double>& Statistics::get_density() {return m_density;}

	/**
	* @brief Accessor for x velocity data
	* @return Vector4D<double>
	* @sa calc_vels()
	*
	* calc_vels() is used to turn the data in this into an actual velocity.
	* Before doing that it's just a collection of Monte Carlo type values.
	*/
	Vectors::Vector4D<double>& Statistics::get_vx() {return m_vx;}

	/**
	* @brief Accessor for y velocity data
	* @return Vector4D<double>
	* @sa calc_vels()
	*
	* calc_vels() is used to turn the data in this into an actual velocity.
	* Before doing that it's just a collection of Monte Carlo type values.
	*/
	Vectors::Vector4D<double>& Statistics::get_vy() {return m_vy;}

	/**
	* @brief Accessor for z velocity data
	* @return Vector4D<double>
	* @sa calc_vels()
	*
	* calc_vels() is used to turn the data in this into an actual velocity.
	* Before doing that it's just a collection of Monte Carlo type values.
	*/
	Vectors::Vector4D<double>& Statistics::get_vz() {return m_vz;}

	/**
	* @brief Accessor for if velocity stats are being tracked
	* @return Returns boolean true if velocity stats are being tracked, false
	* if not.
	*/
	bool Statistics::get_vel_stats() {return m_vel_stats;}

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
		ret_stats.m_vx = m_vx + other.m_vx;
		ret_stats.m_vy = m_vy + other.m_vy;
		ret_stats.m_vz = m_vz + other.m_vz;

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
		const int yidx, const int zidx, const int value)
	{
		m_weights(tidx, xidx, yidx, zidx) += value;
	}

	/**
	*
	*/
	void Statistics::add_vels(const int tidx, const int xidx, const int yidx,
		const int zidx, const double vx, const double vy, const double vz)
	{
		// Add velocities to the running total at each cell location
		m_vx(tidx, xidx, yidx, zidx) += vx;
		m_vy(tidx, xidx, yidx, zidx) += vy;
		m_vz(tidx, xidx, yidx, zidx) += vz;
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

			// Normalize each weight value by the volume and total number of
			// particles launched. This goes from units of (s) to (s/m3). 
			m_density(i,j,k,l) = m_weights(i,j,k,l) / (dx * dy * dz) 
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
		if (!m_vel_stats)
		{
			std::cerr << "Error! Cannot calculate average impurity velocities"
				<< " if velocity tracking is off. Set imp_vel_stats to yes if"
				<< " velocity tracking is desired.\n";
		}
		else
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
				// Convert each velocity sum in each cell to an average velocity
				// by dividing it by the number of counts.
				m_vx(i,j,k,l) /= m_counts(i,j,k,l);
				m_vy(i,j,k,l) /= m_counts(i,j,k,l);
				m_vz(i,j,k,l) /= m_counts(i,j,k,l);
			}
			}
			}	
			}
		}
	}
}
