/**
* @file sputtering.h
*
* @brief Header file for sputtering.cpp
*/

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <cstdint>
#include <string>
#include <unordered_map>
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
	* One piece of the wall surface that the plasma can sputter, in one
	* toroidal cell.
	*
	* On a SOLEDGE3X background it is a face of a plasma cell (chi = 0)
	* whose neighbor across that face is a wall cell (chi = 1). Its
	* position is that face's, but its normal and area come from the
	* nearest segment of the wall outline, unless there isn't one close
	* enough (fallback). On a grid fitted to the targets (HERMES) it is a
	* cell face on a wall face of its zone, with the face's own normal and
	* area, and there is no wall cell (wall_zone = -1). See sputtering.cpp
	* for why.
	*
	* On a Cartesian zone (the straight slab, Zones::Zones::is_cartesian)
	* the element is a cell face on a wall face of its zone, as for HERMES,
	* but with no toroidal angle: the fields named R, Z and phi hold X, Y
	* and Z, so the in-plane geometry is the same code (cartesian = true).
	* Use normal_xyz for the unit normal in Cartesian coordinates.
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

		int wall_zone {};		///< Zone of the wall cell across the face,
								///< or -1 if the face is a wall face of
								///< its zone (a target), with no wall cell
		int wall_i {};			///< x index of the wall cell
		int wall_j {};			///< y index of the wall cell (its z index
								///< is k, the same as the plasma cell's)

		bool cartesian {};		///< On a Cartesian zone: R, Z and phi
								///< above are X, Y and Z

		/// True if the element is on a wall face of its zone (a target)
		bool on_wall_face() const {return wall_zone < 0;}

		/**
		* The unit normal into the plasma in Cartesian (X, Y, Z).
		*
		* On a toroidal zone the normal (n_R, n_Z) turns with the toroidal
		* angle: positions are (R cos phi, phi_sign R sin phi, Z), so R
		* points along (cos phi, phi_sign sin phi, 0). On a Cartesian zone
		* it's (n_R, n_Z, 0) everywhere, since there n_R and n_Z are its X
		* and Y parts and wall faces are x or y faces.
		*
		* @param z The element's z coordinate where the normal is wanted
		*        (the toroidal angle; not used on a Cartesian zone).
		* @param phi_sign The background's phi_sign (+1 or -1).
		* @return The unit normal (X, Y, Z).
		*/
		std::array<double, 3> normal_xyz(const double z,
			const double phi_sign) const
		{
			if (cartesian) return {n_R, n_Z, 0.0};
			return {n_R * std::cos(z), n_R * phi_sign * std::sin(z), n_Z};
		}
	};

	/**
	* The wall surfaces of a background that the plasma can sputter, as a
	* list of SurfaceElements: the staircase between plasma and wall cells
	* of a SOLEDGE3X background, or the target faces of a grid fitted to
	* them (HERMES). See sputtering.cpp for the member function docs.
	*/
	class Surface
	{
	private:
		std::vector<SurfaceElement> m_elements {};	///< The elements
		int m_num_fallback {};		///< Elements that fell back
		double m_total_area {};		///< Sum of the elements' areas (m^2)
		bool m_wall_faces {};		///< Elements are zone wall faces

		void find_staircase(const Background::Background& bkg,
			double max_wall_dist);
		void find_wall_faces(const Background::Background& bkg);

	public:

		Surface() = default;

		explicit Surface(const Background::Background& bkg,
			double max_wall_dist = 2.0);

		const std::vector<SurfaceElement>& elements() const;
		std::size_t size() const;
		int num_fallback() const;
		double total_area() const;
		bool from_wall_faces() const;

		void print_summary(std::ostream& os) const;
	};

	bool has_surfaces(const Background::Background& bkg);

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
		double m_phi_sign {1.0};		///< Background's toroidal handedness
		std::vector<double> m_e_sheath {};	///< Sheath energy (eV), per idx
		std::vector<double> m_ti {};		///< Ion temperature (eV), per idx
		std::vector<double> m_phi_wall {};	///< Floating drop (V), per idx
		std::vector<std::array<double, 4>> m_cell_bounds {};	///< x and y
										///< edges of each element's cell
		std::vector<double> m_frame_lo {};	///< Start of each frame's window
		std::vector<double> m_frame_hi {};	///< End of each frame's window
		std::vector<double> m_cdf {};		///< Cumulative atoms per idx
		double m_source_rate {};			///< Time-averaged total rate
		std::vector<double> m_mean_rate {};	///< Time-averaged rate per element

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
		const std::vector<double>& mean_rate() const;
		int num_above_table() const;
		const std::vector<double>& frame_lo() const;
		const std::vector<double>& frame_hi() const;
		const YieldTable& table() const;
		const std::vector<std::array<double, 4>>& cell_bounds() const;
		const std::vector<double>& phi_wall() const;
		double phi_sign() const;

		Birth sample_birth(pcg32& rng) const;

		void print_summary(std::ostream& os) const;
	};

	/**
	* Where the sputtered atoms end up back on the wall, and how long and
	* in what charge state they took to get there. See sputtering.cpp for
	* the member function docs.
	*
	* Everything is a rate, in atoms/s, since every particle's weight is
	* the atoms/s it stands for: adding up the weights of the particles
	* that land on an element gives the rate atoms are redeposited on it.
	*
	* Prompt returns are also tallied on their own: those landing within
	* one gyroperiod of the singly charged ion after the atom first
	* ionized (Slots::Slots::prompt_by), as in Fussmann's model of prompt
	* redeposition.
	*/
	/**
	* What an absorbed particle hits the wall with: its Cartesian velocity
	* and mass, for the impact energy and angle tallies of Deposition.
	*/
	struct Impact
	{
		std::array<double, 3> v {};	///< Velocity (vX, vY, vZ) (m/s)
		double mass {};				///< Mass (kg)
	};

	class SelfSource;

	class Deposition
	{
	private:
		const Surface* m_surface {nullptr};	///< The surface elements
		SelfSource* m_self {nullptr};	///< Self-sputtering, if on
		std::unordered_map<std::uint64_t, std::vector<std::size_t>>
			m_wall_cells {};		///< Elements bordering each wall cell
		std::unordered_map<std::uint64_t, std::size_t>
			m_wall_faces {};		///< Element on each cell of a wall face
		std::vector<double> m_redeposited {};	///< Rate onto each element
		double m_unassigned {};		///< Into wall cells with no element
		double m_wall_face {};		///< At wall faces with no element
		std::vector<double> m_lifetime_edges {};	///< Time since birth
													///< bin edges (s)
		std::vector<double> m_lifetime {};	///< Rate in each time bin
		std::vector<double> m_charge {};	///< Rate at each charge state
		std::vector<double> m_prompt {};	///< Prompt rate onto each
											///< element
		std::vector<double> m_prompt_charge {};	///< Prompt rate at each
												///< charge state
		double m_phi_sign {1.0};	///< Handedness of (R, phi, Z)
		double m_lifetime_sum {};	///< Sum of weight * time since birth
		std::vector<double> m_energy_edges {};	///< Impact energy bin
												///< edges (eV)
		std::vector<double> m_energy {};	///< Rate in each energy bin
		std::vector<double> m_angle {};		///< Rate in each 1 degree bin
											///< of the angle from the
											///< normal
		double m_energy_sum {};		///< Sum of weight * impact energy
		double m_energy_weight {};	///< Rate with an impact energy
		double m_angle_sum {};		///< Sum of weight * impact angle
		double m_angle_weight {};	///< Rate with an impact angle
		std::vector<double> m_element_energy {};	///< Sum of weight *
													///< impact energy on
													///< each element

		bool add_lifetime_and_charge(double pt, int q, double weight,
			double prompt_by);
		void add_impact(const Impact* impact, const SurfaceElement* el,
			double z, double weight, std::size_t e);

	public:

		Deposition() = default;

		Deposition(const Surface& surface, int atomic_number,
			double phi_sign = 1.0);

		long add_into_wall_cell(int zone, int i, int j, int k,
			int prev_zone, int prev_i, int prev_j, double pt, int q,
			double weight, double prompt_by = -1.0,
			const Impact* impact = nullptr);
		long add_at_wall_face(const Zones::Zones& zones, int zone,
			Zones::Face face, const std::array<double, 3>& pos, double pt,
			int q, double weight, double prompt_by = -1.0,
			const Impact* impact = nullptr);

		void set_self(SelfSource* self);
		SelfSource* self() const;

		void reduce(int root);

		const std::vector<double>& redeposited() const;
		double unassigned() const;
		double wall_face() const;
		double total() const;
		const std::vector<double>& lifetime_edges() const;
		const std::vector<double>& lifetime() const;
		const std::vector<double>& charge() const;
		const std::vector<double>& prompt() const;
		const std::vector<double>& prompt_charge() const;
		double prompt_total() const;
		double mean_lifetime() const;
		const std::vector<double>& energy_edges() const;
		const std::vector<double>& impact_energy() const;
		const std::vector<double>& impact_angle() const;
		double mean_energy() const;
		double mean_angle() const;
		const std::vector<double>& element_energy() const;

		void print_summary(std::ostream& os, double source_rate) const;
	};
	/**
	* The energy and angle an impurity hits the wall with, for the
	* self-sputtering yield (see self_impact).
	*/
	struct SelfImpact
	{
		double energy {};	///< Impact energy (eV)
		double angle {};	///< Angle from the element's normal (degrees)
	};

	SelfImpact self_impact(const std::array<double, 3>& v, double mass,
		const std::array<double, 3>& n, int q, double phi_wall);

	double roulette(double weight, double min_weight, double u);

	/**
	* A particle that landed on the wall and carries on as one of the atoms
	* it sputtered (SelfSource::rebirth). If alive is false it ends there,
	* redeposited, and nothing else is set.
	*/
	struct Rebirth
	{
		bool alive {};			///< Carries on as a sputtered atom
		bool avalanche {};		///< Ended at sput_self_max_gen
		int zone {};			///< Zone it starts in
		double x {};			///< Computational x
		double y {};			///< Computational y
		double z {};			///< Computational z
		int i {};				///< x index of its cell
		int j {};				///< y index of its cell
		int k {};				///< z index of its cell
		double vX {};			///< Cartesian velocity (m/s)
		double vY {};
		double vZ {};
		double weight {};		///< New weight (atoms/s)
	};

	/**
	* Self-sputtering: impurities landing on a wall element sputter more
	* atoms, from a yield table of the impurity on the wall material (e.g.,
	* W on W). The particle carries on as one of those atoms, with its
	* weight multiplied by the yield, so no particles are added mid-run.
	* See sputtering.cpp for how and the member function docs.
	*
	* The tallies are rates (atoms/s), like Deposition's. Per-element
	* values are averaged over the run, not kept per frame.
	*/
	class SelfSource
	{
	private:
		YieldTable m_table {};			///< Impurity on the wall material
		const Source* m_source {nullptr};	///< Elements, cells, drops
		bool m_sheath_on {};			///< Push includes the sheath field
		double m_lambda {};				///< Largest energy fraction the
										///< impurity can give a wall atom
		double m_imp_mass_kg {};		///< Mass of the impurity (kg)
		double m_phi_sign {1.0};		///< Background's toroidal handedness
		double m_min_weight {};			///< Roulette below this (atoms/s)
		int m_max_gen {};				///< Generation that stops the run

		std::vector<double> m_self {};	///< Self-sputtered rate per element
		std::vector<double> m_gen_weight {};	///< Weight carried on into
												///< each generation
		std::vector<double> m_gen_landed {};	///< Weight landing on
												///< elements from each
												///< generation
		double m_landed {};				///< Rate landing on elements
		double m_self_total {};			///< Self-sputtered rate, sum w Y
		double m_energy_sum {};			///< Sum of w E, landed on elements
		double m_angle_sum {};			///< Sum of w angle
		double m_roulette_lost {};		///< Rate ended by roulette
		long long m_num_below {};		///< Landings below threshold
		long long m_num_above {};		///< Landings past the table's end

	public:

		SelfSource() = default;

		SelfSource(const Source& source, const YieldTable& table,
			const Options::Options& opts, bool sheath_on);

		Rebirth rebirth(std::size_t e, int frame,
			const std::array<double, 3>& pos, const Impact& impact, int q,
			double weight, int generation, pcg32& rng);

		void reduce(int root);

		const YieldTable& table() const;
		double min_weight() const;
		int max_gen() const;
		const std::vector<double>& self_sputtered() const;
		const std::vector<double>& generation_weight() const;
		const std::vector<double>& generation_landed() const;
		double self_total() const;
		double landed() const;
		double mean_yield() const;
		double mean_energy() const;
		double mean_angle() const;
		double roulette_lost() const;
		double chain_multiplication() const;
		long long num_below_threshold() const;
		long long num_above_table() const;

		void print_summary(std::ostream& os, double source_rate) const;
	};
}
