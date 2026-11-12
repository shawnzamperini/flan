/**
* @file sputtering.h
*
* @brief Header file for sputtering.cpp
*/

#pragma once

#include <string>
#include <vector>


namespace Sputtering
{
	/**
	* A sputtering yield table, Y(E, angle), made by
	* python/rustbca_yields.py. See sputtering.cpp for the member function
	* docs.
	*/
	class YieldTable
	{
	private:
		std::string m_path {};				///< File the table was read from
		std::vector<double> m_log_energy {};	///< log10(energy [eV])
		std::vector<double> m_angle {};		///< Angle from the normal (deg)
		std::vector<double> m_yield {};		///< Yield, (energy, angle) row major
		double m_energy_min {};			///< Lowest energy (eV)
		double m_energy_max {};			///< Highest energy (eV)
		double m_projectile_Z {};		///< Projectile atomic number
		double m_projectile_mass_amu {};	///< Projectile mass (amu)
		double m_target_Z {};			///< Target atomic number
		double m_target_mass_amu {};	///< Target mass (amu)
		double m_target_Es {};			///< Target surface binding energy (eV)

	public:

		YieldTable() = default;

		explicit YieldTable(const std::string& path);

		double get_yield(double energy, double angle) const;

		const std::string& path() const;
		double energy_min() const;
		double energy_max() const;
		double projectile_Z() const;
		double projectile_mass_amu() const;
		double target_Z() const;
		double target_mass_amu() const;
		double target_Es() const;
	};
}
