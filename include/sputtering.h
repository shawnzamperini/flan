/**
* @file sputtering.h
*
* @brief Header file for sputtering.cpp
*/

#pragma once

#include <array>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "background.h"
#include "options.h"
#include "pcg32.h"
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

	double sheath_energy(double te, double ti, double elec_mass_amu,
		double ion_mass_amu, int charge);

	/// An effective yield and how much of the energy distribution was past
	/// the end of the yield table
	struct EffectiveYield
	{
		double yield {};		///< Yield averaged over the energies
		double frac_above {};	///< Fraction of ions above the table
	};

	EffectiveYield effective_yield(const YieldTable& table, double e_sheath,
		double ti, double angle);

	double sample_impact_energy(const YieldTable& table, double e_sheath,
		double ti, double angle, double u);

	double sample_thompson(double e_b, double e_max, double u);

	std::array<double, 3> cosine_direction(const std::array<double, 3>& n,
		double u1, double u2);

	/**
	* Where, when and how fast a sputtered atom starts (see
	* Source::sample_birth), plus which frame and element it came from and
	* the energy of the ion that sputtered it, for checking.
	*/
	struct Birth
	{
		int zone {};			///< Zone it starts in
		double t {};			///< Time (s)
		double x {};			///< Computational x
		double y {};			///< Computational y
		double z {};			///< Computational z (phi for SOLEDGE3X)
		double vX {};			///< Cartesian velocity (m/s)
		double vY {};
		double vZ {};
		int frame {};			///< Time frame it was picked from
		std::size_t element {};	///< Surface element it was picked from
		double e_impact {};		///< Energy of the sputtering ion (eV)
		double e_atom {};		///< Energy of the sputtered atom (eV)
	};

	/**
	* The rate atoms are sputtered off each surface element by the
	* background plasma, in each background time frame. See sputtering.cpp
	* for how it's calculated and the member function docs.
	*
	* Per-element values are stored for every frame, frame by frame: the
	* value for frame t and element e is at index t * size() + e of the
	* vectors.
	*/
	class Source
	{
	private:
		Surface m_surface {};				///< The surface elements
		int m_nframes {};					///< Number of time frames
		std::vector<double> m_flux {};		///< Ion flux onto the wall
		std::vector<double> m_flux_par {};	///< Parallel part of the flux
		std::vector<double> m_flux_exb {};	///< E x B part of the flux
		std::vector<double> m_energy {};	///< Mean impact energy (eV)
		std::vector<double> m_yield {};		///< Effective yield
		std::vector<double> m_rate {};		///< Sputtered atoms/s
		std::vector<double> m_total_rate {};	///< Sum over elements, per frame
		int m_num_above_table {};	///< Element-frames past the table
		std::vector<std::string> m_warnings {};	///< For print_summary

		// For starting particles (see sample_birth)
		YieldTable m_table {};				///< The yield table
		double m_angle {};					///< Ion angle of incidence (deg)
		double m_lambda {};				///< Largest energy fraction an ion
										///< can give a target atom
		double m_imp_mass_kg {};		///< Mass of the sputtered atoms
		std::vector<double> m_e_sheath {};	///< Sheath energy (eV), per idx
		std::vector<double> m_ti {};		///< Ion temperature (eV), per idx
		std::vector<std::array<double, 4>> m_cell_bounds {};	///< x and y
										///< edges of each element's cell
		std::vector<double> m_frame_lo {};	///< Start of each frame's window
		std::vector<double> m_frame_hi {};	///< End of each frame's window
		std::vector<double> m_cdf {};		///< Cumulative atoms per idx
		double m_source_rate {};			///< Time-averaged total rate

	public:

		Source() = default;

		Source(const Surface& surface, const Background::Background& bkg,
			const YieldTable& table, const Options::Options& opts);

		const Surface& surface() const;
		std::size_t size() const;
		int nframes() const;
		const std::vector<double>& flux() const;
		const std::vector<double>& flux_par() const;
		const std::vector<double>& flux_exb() const;
		const std::vector<double>& energy() const;
		const std::vector<double>& yield() const;
		const std::vector<double>& rate() const;
		const std::vector<double>& total_rate() const;
		double mean_total_rate() const;
		double source_rate() const;
		int num_above_table() const;
		const std::vector<double>& frame_lo() const;
		const std::vector<double>& frame_hi() const;

		Birth sample_birth(pcg32& rng) const;

		void print_summary(std::ostream& os) const;
	};
}
