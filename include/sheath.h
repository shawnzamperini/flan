#pragma once

/**
* @file  sheath.h
*
* @brief Header file for sheath.cpp: the electric field of the sheath in
*        front of a wall
*/
#include <array>
#include <cstdint>
#include <iosfwd>
#include <unordered_map>
#include <vector>

#include "background.h"
#include "sputtering.h"

namespace Sheath
{
	/// Which profile the potential has (the sheath_model option)
	enum class Model : int
	{
		off = 0,		///< No sheath field
		brooks = 1,		///< Brooks' two exponentials
		borodkina = 2	///< Borodkina and Tsvetkov's analytic profile
	};

	/**
	* The sheath in front of one wall element in one time frame: what the
	* potential and its field at a distance from the wall depend on. Made
	* by make_params from the plasma next to the wall.
	*
	* Potentials are relative to the plasma at the entrance of the
	* magnetic presheath, so they are 0 far from the wall and phi_wall
	* (negative) at it. lam_* are those potentials over Te (e phi / Te),
	* and xi_* and L_mps distances over the Debye length, as in Borodkina
	* and Tsvetkov (2015).
	*/
	struct Params
	{
		double te {};		///< Electron temperature (eV)
		double phi_wall {};	///< Potential at the wall, the floating drop (V)
		double f_d {};		///< Share of the drop in the Debye sheath
		double lambda_d {};	///< Electron Debye length (m)
		double rho_s {};	///< Ion sound gyroradius c_s / Omega_i (m)
		double theta {};	///< Angle between B and the wall's normal (rad)

		// Borodkina and Tsvetkov's profile
		double lam_w {};	///< e phi_wall / Te
		double lam_mps {};	///< e phi / Te at the presheath/Debye boundary
		double a {};		///< Decay rate of the Debye layer (1 / lambda_d)
		double q {};		///< Amplitude of the Debye layer (Q)
		double xi_mps {};	///< Presheath/Debye boundary (lambda_d)
		double l_mps {};	///< Decay length of the presheath (lambda_d)
		bool debye {};		///< There is a Debye layer (below the critical
							///< angle)
	};

	Params make_params(double te, double ti, double ne, double bmag,
		double cos_theta, double elec_mass_amu, double ion_mass_amu);

	double potential(const Params& p, Model model, double d);

	double field(const Params& p, Model model, double d);

	double range(const Params& p, Model model, double n_scales = 28.0);

	/**
	* The sheath field in front of the wall surfaces of a background: the
	* Params of every surface element in every time frame, and for each
	* plasma cell next to the wall, its elements. field_at gives the field
	* at a particle. See sheath.cpp for the member function docs.
	*/
	class Sheath
	{
	public:
		/// Most a sub-step may move a particle, in local field scale
		/// lengths
		static constexpr double substep_frac {0.02};

		/// Most sub-steps in one time step
		static constexpr int max_substeps {20000};

	private:
		Model m_model {Model::off};
		int m_nframes {};
		double m_phi_sign {1.0};

		/// Per element: the plasma cell's zone, face and the cell's extent
		/// along the face's axis, the cell's thickness along the normal
		/// (m), the normal (n_R, n_Z or n_X, n_Y) and whether it's
		/// Cartesian
		struct Element
		{
			int zone {};
			Zones::Face face {};
			double lo {};
			double hi {};
			double thickness {};
			double n_R {};
			double n_Z {};
			bool cartesian {};
		};
		std::vector<Element> m_elements {};

		/// Params for each frame and element, frame-major
		std::vector<Params> m_params {};

		/// Elements of each plasma cell next to the wall
		std::unordered_map<std::uint64_t, std::vector<std::size_t>>
			m_cells {};

		int m_num_thin {};	///< Elements whose cell is thinner than the
							///< sheath
		double m_max_lost {};	///< Most of the drop beyond a cell's far
								///< side, as a fraction

	public:
		Sheath() = default;

		Sheath(const Sputtering::Surface& surface,
			const Background::Background& bkg, Model model);

		Model model() const;
		std::size_t size() const;
		int nframes() const;
		const Params& params(int frame, std::size_t element) const;
		int num_thin() const;
		double max_lost() const;

		bool distance(int zone, int i, int j, int k,
			const std::array<double, 3>& pos, std::size_t& element,
			double& d) const;

		std::array<double, 3> field_at(int zone, int frame, int i, int j,
			int k, const std::array<double, 3>& pos) const;

		int substeps(int zone, int frame, int i, int j, int k,
			const std::array<double, 3>& pos, double speed, double q_over_m,
			double dt) const;

		void print_summary(std::ostream& os) const;
	};
}
