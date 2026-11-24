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
		double ti {};		///< Ion temperature (eV)
		double bmag {};		///< Field strength (T)
		double ion_mass_amu {};	///< Main ion mass (amu)
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
		double cos_theta, double elec_mass_amu, double ion_mass_amu,
		double f_d = -1.0, double drop_te = -1.0);

	double potential(const Params& p, Model model, double d);

	double field(const Params& p, Model model, double d);

	double range(const Params& p, Model model, double n_scales = 28.0);

	int substeps(const Params& p, Model model, double d, double speed,
		double q_over_m, double dt);

	/// One test ion that hit the wall (sample_impacts)
	struct ImpactSample
	{
		double energy {};	///< Kinetic energy when it hit (eV)
		double angle {};	///< Angle of its velocity from the normal
							///< (degrees, 0 to 90)
		double weight {};	///< Its share of the flux into the sheath:
							///< its speed along B toward the wall at the
							///< start (m/s)
	};

	/// What sample_impacts returns
	struct Impacts
	{
		std::vector<ImpactSample> hits {};	///< The ions that hit the wall
		int num_ions {};			///< Ions started
		int lost {};				///< Ions that went back out into the
									///< plasma or ran out of steps
		double start_distance {};	///< Where they started (m)
		double max_energy_error {};	///< Most |E - E_start - (phi_start -
									///< phi_wall)| over the hits, over
									///< the drop
	};

	Impacts sample_impacts(const Params& p, Model model, int num_ions,
		std::uint64_t seed, double start_rho_s = 10.0,
		int steps_per_gyro = 100, int max_gyros = 2000);

	/**
	* D+ impact samples (sample_impacts) on a grid of the three numbers
	* the one-dimensional problem depends on, for the sputtered source
	* (sput_impact_opt = "sheath"). In units of rho_s, c_s and Te, an
	* ion's path through the sheath depends only on the field angle theta
	* (which also sets f_D), Ti / Te and lambda_D / rho_s, so impact
	* energies scale with Te and the angles don't change. The grid is
	* theta_step in theta, a factor exp(log_ti_step) in Ti / Te and
	* 2^log2_lambda_step in lambda_D / rho_s. Only the points the run
	* uses are made (add), then sampled (build). See sheath.cpp.
	*/
	class ImpactTable
	{
	public:
		static constexpr double theta_step {0.5};		///< Degrees
		static constexpr double log_ti_step {0.1};		///< ln(Ti / Te)
		static constexpr double log2_lambda_step {1.0};	///< log2(lambda_D
														///< / rho_s)
		static constexpr int default_num_ions {2000};	///< Per point

	private:
		Model m_model {Model::off};
		double m_me_amu {};
		double m_mi_amu {};
		double m_f_d {-1.0};
		double m_drop_te {-1.0};
		int m_num_ions {default_num_ions};
		std::uint64_t m_seed {1};
		std::unordered_map<std::uint64_t, int> m_index {};	///< Key to point
		std::vector<std::array<int, 3>> m_keys {};	///< Grid indices of
													///< each point
		/// Each point's samples, energies over Te and weights summing to 1
		std::vector<std::vector<ImpactSample>> m_samples {};
		int m_lost {};
		double m_max_energy_error {};
		bool m_built {};

	public:
		ImpactTable() = default;

		ImpactTable(Model model, double elec_mass_amu, double ion_mass_amu,
			double f_d = -1.0, double drop_te = -1.0,
			int num_ions = default_num_ions, std::uint64_t seed = 1);

		int add(const Params& p);
		void build();

		Model model() const;
		std::size_t size() const;
		int num_ions() const;
		bool built() const;
		const std::vector<ImpactSample>& samples(int point) const;
		std::array<double, 3> point_values(int point) const;
		Params point_params(int point) const;
		int lost() const;
		double max_energy_error() const;
	};

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
			const Background::Background& bkg, Model model,
			double f_d = -1.0, double drop_te = -1.0);

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
