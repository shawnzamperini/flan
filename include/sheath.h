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
#include "sheath_hd.h"
#include "sputtering.h"

namespace Sheath
{

	Params make_params(double te, double ti, double ne, double bmag,
		double cos_theta, double elec_mass_amu, double ion_mass_amu,
		double f_d = -1.0, double drop_te = -1.0);

	// Model, Params, potential, field, range and substeps are in
	// sheath_hd.h, shared with the GPU

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
		static constexpr double substep_frac {::Sheath::substep_frac};

		/// Most sub-steps in one time step
		static constexpr int max_substeps {::Sheath::max_substeps};

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

		/// The elements as flat arrays, for SheathView
		struct Flat
		{
			std::vector<int> zone {};
			std::vector<int> face {};
			std::vector<double> lo {};
			std::vector<double> hi {};
			std::vector<double> thickness {};
			std::vector<double> n_R {};
			std::vector<double> n_Z {};
			std::vector<int> cartesian {};
		} m_flat {};

		/// Plasma cells next to the wall, sorted by cell_key, and the
		/// elements of each (m_cell_start[c] to [c + 1] in m_cell_elems,
		/// in element order). Sorted arrays rather than a hash map, so the
		/// GPU searches the same ones.
		std::vector<std::uint64_t> m_cell_keys {};
		std::vector<int> m_cell_start {};
		std::vector<int> m_cell_elems {};

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

		SheathView view() const;

		SheathDevice to_device(int device_id) const;
	};

	void free_sheath(SheathDevice& sheath_d);
}
