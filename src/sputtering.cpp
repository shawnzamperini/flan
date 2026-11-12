/**
* @file sputtering.cpp
*
* @brief Sputtering of wall surfaces by the background plasma.
*
* The sputtering yield Y (sputtered atoms per incident ion) depends on the
* incident ion's energy E and angle from the surface normal. It is read
* from a table made with RustBCA by python/rustbca_yields.py (see
* data/sputtering/) and bilinearly interpolated, linear in log10(E) and in
* angle.
*/

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "errors.h"
#include "hdf5_utils.h"
#include "sputtering.h"


namespace Sputtering
{
	namespace
	{
		/**
		* Finds the interval of a sorted grid to interpolate in.
		*
		* @param grid Strictly increasing values, at least two.
		* @param value Value to find. Must be within the grid.
		* @return Index i of the lower end of the interval, so the value is
		*         between grid[i] and grid[i + 1], from 0 to size - 2.
		*/
		std::size_t lower_index(const std::vector<double>& grid,
			const double value)
		{
			const auto upper {std::upper_bound(grid.begin(), grid.end(),
				value)};
			const std::size_t i {static_cast<std::size_t>(
				std::max<std::ptrdiff_t>(upper - grid.begin() - 1, 0))};
			return std::min(i, grid.size() - 2);
		}

		/**
		* Checks that a grid has at least two values and is strictly
		* increasing, stopping the run if not.
		*
		* @param grid Values to check.
		* @param name Name of the grid, for the error message.
		* @param path File the grid is from, for the error message.
		*/
		void check_grid(const std::vector<double>& grid,
			const std::string& name, const std::string& path)
		{
			if (grid.size() < 2)
			{
				Errors::fatal_error("The sputtering yield table " + path
					+ " has " + std::to_string(grid.size()) + " " + name
					+ " value(s). At least 2 are needed to interpolate.");
			}
			for (std::size_t i {1}; i < grid.size(); ++i)
			{
				if (!(grid[i] > grid[i - 1]))
				{
					Errors::fatal_error("The " + name + " values in the "
						"sputtering yield table " + path + " are not "
						"strictly increasing.");
				}
			}
		}
	}

	/**
	* Reads a sputtering yield table.
	*
	* The file is an HDF5 file made by python/rustbca_yields.py, with
	* /energy (eV), /angle (degrees from the surface normal), /yield (atoms
	* per ion, shape (energy, angle)), /projectile/Z, /projectile/mass_amu,
	* /target/Z, /target/mass_amu and /target/Es (eV). Its other datasets
	* and attributes aren't needed. The run stops if a dataset is missing,
	* the yield has the wrong shape, there are fewer than two energies or
	* angles, the energies aren't positive and strictly increasing, the
	* angles aren't strictly increasing between 0 and 90, or a yield is
	* negative or not finite.
	*
	* @param path Path to the table.
	*/
	YieldTable::YieldTable(const std::string& path)
		: m_path {path}
	{
		const Hdf5::File file(path);
		const std::string where {"The sputtering yield table " + path};

		const std::vector<double> energy {Hdf5::read_dataset<double>(file,
			"/energy").data};
		m_angle = Hdf5::read_dataset<double>(file, "/angle").data;
		const Hdf5::Dataset<double> yield {Hdf5::read_dataset<double>(file,
			"/yield")};

		check_grid(energy, "energy", path);
		check_grid(m_angle, "angle", path);
		if (energy.front() <= 0.0)
		{
			Errors::fatal_error(where + " has energies <= 0.");
		}
		if (m_angle.front() < 0.0 || m_angle.back() > 90.0)
		{
			Errors::fatal_error(where + " has angles outside 0 to 90 "
				"degrees.");
		}

		const std::vector<hsize_t> shape {energy.size(), m_angle.size()};
		if (yield.shape != shape)
		{
			Errors::fatal_error(where + " has a yield with the wrong shape. "
				"It should be (" + std::to_string(shape[0]) + ", "
				+ std::to_string(shape[1]) + "), the number of energies "
				"and angles.");
		}
		for (const double y : yield.data)
		{
			if (!std::isfinite(y) || y < 0.0)
			{
				Errors::fatal_error(where + " has a negative or non-finite "
					"yield.");
			}
		}
		m_yield = yield.data;

		m_log_energy.reserve(energy.size());
		for (const double e : energy) m_log_energy.push_back(std::log10(e));
		m_energy_min = energy.front();
		m_energy_max = energy.back();

		m_projectile_Z = Hdf5::read_scalar<double>(file, "/projectile/Z");
		m_projectile_mass_amu = Hdf5::read_scalar<double>(file,
			"/projectile/mass_amu");
		m_target_Z = Hdf5::read_scalar<double>(file, "/target/Z");
		m_target_mass_amu = Hdf5::read_scalar<double>(file,
			"/target/mass_amu");
		m_target_Es = Hdf5::read_scalar<double>(file, "/target/Es");
	}

	/**
	* Returns the sputtering yield at an incident energy and angle.
	*
	* The yield is bilinearly interpolated, linear in log10(E) and in
	* angle. Y itself is interpolated, not log(Y), since the yield is 0
	* below threshold. Below the lowest energy in the table the yield is 0,
	* since the table starts below threshold. Above the highest energy, the
	* highest energy is used. The angle is clamped to the table's angles
	* (e.g., 90 uses the largest angle, 89).
	*
	* @param energy Incident energy (eV).
	* @param angle Incident angle from the surface normal (degrees).
	* @return Sputtering yield (atoms per ion).
	*/
	double YieldTable::get_yield(double energy, double angle) const
	{
		if (!(energy >= m_energy_min)) return 0.0;

		const double log_energy {std::min(std::log10(energy),
			m_log_energy.back())};
		angle = std::clamp(angle, m_angle.front(), m_angle.back());

		const std::size_t i {lower_index(m_log_energy, log_energy)};
		const std::size_t j {lower_index(m_angle, angle)};
		const double fe {(log_energy - m_log_energy[i])
			/ (m_log_energy[i + 1] - m_log_energy[i])};
		const double fa {(angle - m_angle[j])
			/ (m_angle[j + 1] - m_angle[j])};

		const std::size_t na {m_angle.size()};
		const double y00 {m_yield[i * na + j]};
		const double y01 {m_yield[i * na + j + 1]};
		const double y10 {m_yield[(i + 1) * na + j]};
		const double y11 {m_yield[(i + 1) * na + j + 1]};

		return (1.0 - fe) * ((1.0 - fa) * y00 + fa * y01)
			+ fe * ((1.0 - fa) * y10 + fa * y11);
	}

	/**
	* @name Accessors
	* The file the table was read from, its energy range (eV), and the
	* projectile and target it's for: atomic numbers, masses (amu) and the
	* target's surface binding energy (eV).
	*/
	///@{
	const std::string& YieldTable::path() const {return m_path;}
	double YieldTable::energy_min() const {return m_energy_min;}
	double YieldTable::energy_max() const {return m_energy_max;}
	double YieldTable::projectile_Z() const {return m_projectile_Z;}
	double YieldTable::projectile_mass_amu() const
		{return m_projectile_mass_amu;}
	double YieldTable::target_Z() const {return m_target_Z;}
	double YieldTable::target_mass_amu() const {return m_target_mass_amu;}
	double YieldTable::target_Es() const {return m_target_Es;}
	///@}
}
