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
#include "vectors.h"

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
		// Issue: This code seems to not be correct, not sure how yet...
		// These aren't always needed, can cut back on memory usage by not
		// including them by default.
		//std::cout << "Allocating velocity arrays...\n";
		m_vX.move_into_data(Vectors::Vector4D<BkgFPType> {dim1, dim2, dim3, 
			dim4});
		m_vY.move_into_data(Vectors::Vector4D<BkgFPType> {dim1, dim2, dim3, 
			dim4});
		m_vZ.move_into_data(Vectors::Vector4D<BkgFPType> {dim1, dim2, dim3, 
			dim4});
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
	* @brief Accessor for x velocity data
	* @return Vector4D<BkgFPType>
	* @sa calc_vels()
	*
	* calc_vels() is used to turn the data in this into an actual velocity.
	* Before doing that it's just a collection of Monte Carlo type values.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_vX() {return m_vX;}

	/**
	* @brief Accessor for y velocity data
	* @return Vector4D<BkgFPType>
	* @sa calc_vels()
	*
	* calc_vels() is used to turn the data in this into an actual velocity.
	* Before doing that it's just a collection of Monte Carlo type values.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_vY() {return m_vY;}

	/**
	* @brief Accessor for z velocity data
	* @return Vector4D<BkgFPType>
	* @sa calc_vels()
	*
	* calc_vels() is used to turn the data in this into an actual velocity.
	* Before doing that it's just a collection of Monte Carlo type values.
	*/
	Vectors::Vector4D<BkgFPType>& Statistics::get_vZ() {return m_vZ;}

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
		ret_stats.m_vX = m_vX + other.m_vX;
		ret_stats.m_vY = m_vY + other.m_vY;
		ret_stats.m_vZ = m_vZ + other.m_vZ;

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
		const BkgFPType vZ, const Background::Background& bkg)
	{
		// Add velocities to the running total at each cell location
		m_vX(tidx, xidx, yidx, zidx) += vX;
		m_vY(tidx, xidx, yidx, zidx) += vY;
		m_vZ(tidx, xidx, yidx, zidx) += vZ;

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
			// Convert each velocity sum in each cell to an average velocity
			// by dividing it by the number of counts.
			double counts {static_cast<double>(m_counts(i,j,k,l))};
			if (counts > 0)
			{
				m_vX(i,j,k,l) /= counts;
				m_vY(i,j,k,l) /= counts;
				m_vZ(i,j,k,l) /= counts;
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
}
