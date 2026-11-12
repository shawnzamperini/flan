/**
* @file sputtering.h
*
* @brief Header file for sputtering.cpp
*/

#pragma once

#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "background.h"
#include "zones.h"


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

	/**
	* One piece of the wall surface that the plasma can sputter.
	*
	* It is a face of a plasma cell (chi = 0) whose neighbor across that
	* face is a wall cell (chi = 1), in one toroidal cell. Its position is
	* that face's, but its normal and area come from the nearest segment of
	* the wall outline, unless there isn't one close enough (fallback). See
	* sputtering.cpp for why.
	*/
	struct SurfaceElement
	{
		int zone {};			///< Zone of the plasma cell
		int i {};				///< x (psi) index of the plasma cell
		int j {};				///< y (theta) index of the plasma cell
		int k {};				///< z (phi) index of the plasma cell
		Zones::Face face {};	///< Face of the plasma cell on the wall

		double R0 {};			///< R of the face's first corner (m)
		double Z0 {};			///< Z of the face's first corner (m)
		double R1 {};			///< R of the face's second corner (m)
		double Z1 {};			///< Z of the face's second corner (m)
		double R_mid {};		///< R of the middle of the face (m)
		double Z_mid {};		///< Z of the middle of the face (m)
		double phi_min {};		///< Toroidal angle at one side (rad)
		double phi_max {};		///< Toroidal angle at the other side (rad)

		double n_R {};			///< R part of the unit normal into the plasma
		double n_Z {};			///< Z part of the unit normal into the plasma
		double area_face {};	///< Area of the cell face itself (m^2)
		double area {};			///< Area of the wall it stands for (m^2)
		double wall_dist {};	///< Distance to the wall outline (m)
		bool fallback {};		///< No wall segment close enough was found
	};

	/**
	* The wall surfaces of a SOLEDGE3X background that the plasma can
	* sputter, as a list of SurfaceElements. See sputtering.cpp for the
	* member function docs.
	*/
	class Surface
	{
	private:
		std::vector<SurfaceElement> m_elements {};	///< The elements
		int m_num_fallback {};		///< Elements that fell back
		double m_total_area {};		///< Sum of the elements' areas (m^2)

	public:

		Surface() = default;

		explicit Surface(const Background::Background& bkg,
			double max_wall_dist = 2.0);

		const std::vector<SurfaceElement>& elements() const;
		std::size_t size() const;
		int num_fallback() const;
		double total_area() const;

		void print_summary(std::ostream& os) const;
	};
}
