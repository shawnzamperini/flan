#pragma once
/**
* @file sheath_hd.h
*
* @brief The sheath field in front of the walls, shared by the CPU and the
*        GPU (see host_device.h).
*
* The profiles (potential, field, range), the number of near-wall
* sub-steps, and the lookup of a particle's wall element and distance
* (distance, field_at, substeps on a SheathView) are written once here.
* Sheath::Sheath (sheath.cpp) calls them on a view of its own data, and the
* GPU kernels on a device copy of the same data (Sheath::to_device), so the
* two push particles through the sheath the same way.
*
* A view holds plain pointers and sizes, never std::vector, so it can be
* copied to a kernel as it is.
*/

#include <cmath>
#include <cstddef>
#include <cstdint>

#include "host_device.h"
#include "zones_hd.h"

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

	/// Most a sub-step may move a particle, in local field scale lengths
	inline constexpr double substep_frac {0.02};

	/// Most sub-steps in one time step
	inline constexpr int max_substeps {20000};

	/**
	* The sheath potential at a distance d from the wall (V), relative to
	* the plasma at the presheath entrance: phi_wall at d = 0 and 0 far
	* away.
	*
	* @param p The sheath (make_params).
	* @param model "brooks" or "borodkina"; "off" gives 0.
	* @param d Distance from the wall along its normal (m), >= 0.
	* @return The potential (V).
	*/
	FLAN_HD inline double potential(const Params& p, const Model model,
		const double d)
	{
		if (model == Model::off || p.te == 0.0) return 0.0;
		if (model == Model::brooks)
		{
			const double debye {p.f_d * std::exp(-d / (2.0 * p.lambda_d))};
			const double mps {p.f_d < 1.0
				? (1.0 - p.f_d) * std::exp(-d / p.rho_s) : 0.0};
			return p.phi_wall * (debye + mps);
		}

		// Borodkina
		const double xi {d / p.lambda_d};
		if (p.debye)
		{
			if (xi <= p.xi_mps)
				return p.te * (p.lam_w + p.q * (1.0 - std::exp(-p.a * xi)));
			if (!(p.l_mps > 0.0)) return 0.0;
			return p.te * p.lam_mps * std::exp(-(xi - p.xi_mps) / p.l_mps);
		}
		if (!(p.l_mps > 0.0)) return 0.0;
		return p.te * p.lam_w * std::exp(-xi / p.l_mps);
	}

	/**
	* The sheath field along the wall's normal into the plasma at a
	* distance d from the wall, E . n = -dphi/dd (V/m). It is negative:
	* the field points into the wall.
	*
	* @param p The sheath (make_params).
	* @param model "brooks" or "borodkina"; "off" gives 0.
	* @param d Distance from the wall along its normal (m), >= 0.
	* @return E . n (V/m).
	*/
	FLAN_HD inline double field(const Params& p, const Model model,
		const double d)
	{
		if (model == Model::off || p.te == 0.0) return 0.0;
		if (model == Model::brooks)
		{
			const double l_ds {2.0 * p.lambda_d};
			const double debye {p.f_d / l_ds * std::exp(-d / l_ds)};
			const double mps {p.f_d < 1.0
				? (1.0 - p.f_d) / p.rho_s * std::exp(-d / p.rho_s) : 0.0};
			return p.phi_wall * (debye + mps);
		}

		// Borodkina: E . n = -(Te / lambda_D) dlam/dxi
		const double xi {d / p.lambda_d};
		double dlam {};
		if (p.debye && xi <= p.xi_mps)
			dlam = p.a * p.q * std::exp(-p.a * xi);
		else if (p.l_mps > 0.0)
		{
			const double lam0 {p.debye ? p.lam_mps : p.lam_w};
			const double xi0 {p.debye ? p.xi_mps : 0.0};
			dlam = -lam0 / p.l_mps * std::exp(-(xi - xi0) / p.l_mps);
		}
		return -p.te * dlam / p.lambda_d;
	}

	/**
	* How far from the wall the field reaches: n_scales of its longest
	* scale length (past the Debye layer, for "borodkina"). With the
	* default 28 the field beyond is below exp(-28) ~ 7e-13 of the drop per
	* scale length, and is left out.
	*
	* @param p The sheath (make_params).
	* @param model The profile.
	* @param n_scales How many scale lengths.
	* @return The distance (m).
	*/
	FLAN_HD inline double range(const Params& p, const Model model,
		const double n_scales = 28.0)
	{
		if (model == Model::off || p.te == 0.0) return 0.0;
		const double l_ds {2.0 * p.lambda_d};
		if (model == Model::brooks)
		{
			return n_scales * (p.f_d < 1.0 ? std::fmax(l_ds, p.rho_s)
				: l_ds);
		}
		// Borodkina: the Debye layer then the presheath, or one of them
		if (p.debye && !std::isfinite(p.xi_mps))
			return n_scales / p.a * p.lambda_d;
		if (p.debye) return (p.xi_mps + n_scales * p.l_mps) * p.lambda_d;
		return n_scales * p.l_mps * p.lambda_d;
	}

	/**
	* Number of sub-steps a particle at distance d from a wall needs this
	* time step, so that none moves it more than substep_frac of the
	* field's local scale length.
	*
	* A particle that can't reach the field this step (or has no charge)
	* needs 1, a normal step; one that can, at least 2. Otherwise the scale
	* length is the Debye layer's (2 lambda_D for "brooks", lambda_D / a
	* for "borodkina") if the particle can reach the Debye layer this step,
	* and the presheath's (rho_s, or rho_s sin(theta)) if not. The speed it
	* might reach is its own plus what the rest of the drop could add,
	* sqrt(2 (q/m) (phi(d) - phi_wall)), so the count holds as it speeds up
	* toward the wall. The count is capped at max_substeps. The main loop
	* (through Sheath::Sheath::substeps or the GPU's SheathView version)
	* and the impact sampler (sample_impacts) both use it, so they
	* sub-step the same way.
	*
	* @param p The sheath (make_params).
	* @param model "brooks" or "borodkina"; "off" gives 1.
	* @param d Distance from the wall along its normal (m).
	* @param speed The particle's speed (m/s).
	* @param q_over_m Its charge over mass (C/kg), > 0 for positive ions.
	* @param dt The time step (s).
	* @return The number of sub-steps, 1 if it needs none.
	*/
	FLAN_HD inline int substeps(const Params& p, const Model model,
		const double d, const double speed, const double q_over_m,
		const double dt)
	{
		if (model == Model::off || !(q_over_m > 0.0) || p.te == 0.0)
			return 1;

		// The fastest it could be going, and how close it could get
		const double drop {std::fmax(potential(p, model, d) - p.phi_wall,
			0.0)};
		const double v_max {speed + std::sqrt(2.0 * q_over_m * drop)};
		const double reach {v_max * dt};
		const double nearest {d - reach};
		if (nearest > range(p, model, 7.0)) return 1;

		// The Debye layer's scale and extent, and the presheath's scale
		double l_debye {};
		double debye_edge {};
		double l_mps {};
		if (model == Model::brooks)
		{
			l_debye = p.f_d > 0.0 ? 2.0 * p.lambda_d : 0.0;
			debye_edge = 7.0 * l_debye;
			l_mps = p.f_d < 1.0 ? p.rho_s : 0.0;
		}
		else
		{
			l_debye = p.debye ? p.lambda_d / p.a : 0.0;
			debye_edge = p.debye ? (std::isfinite(p.xi_mps) ? p.xi_mps
				* p.lambda_d : 0.0) + 7.0 * l_debye : 0.0;
			l_mps = p.l_mps * p.lambda_d;
		}
		double scale {};
		if (l_debye > 0.0 && nearest <= debye_edge) scale = l_debye;
		else scale = l_mps;
		if (l_debye > 0.0 && l_mps > 0.0 && nearest <= debye_edge)
			scale = std::fmin(l_debye, l_mps);
		if (!(scale > 0.0)) return 1;

		// At least 2 in the field, so that a particle crossing the wall
		// gets its impact velocity (see Boris::substep_cpu)
		const double n {std::ceil(reach / (substep_frac * scale))};
		const double capped {std::fmin(std::fmax(n, 2.0),
			static_cast<double>(max_substeps))};
		return static_cast<int>(capped);
	}

	/**
	* A Sheath::Sheath as flat arrays (see that class): per wall element
	* the plasma cell's zone, the face, the cell's extent across the face,
	* its thickness along the element's normal and the normal; the Params
	* of every frame and element (frame-major); and the plasma cells next
	* to the wall, sorted by key, with each one's elements
	* (cell_start[c] to cell_start[c + 1] in cell_elems, in element
	* order).
	*/
	struct SheathView
	{
		Model model {Model::off};		///< The profile
		int nframes {};					///< Time frames
		double phi_sign {1.0};			///< Handedness of (R, phi, Z)
		std::size_t n_elem {};			///< Wall elements
		const int* zone {};				///< Zone of each element's cell
		const int* face {};				///< Zones::Face of each element
		const double* lo {};			///< Cell's lower edge across the face
		const double* hi {};			///< Its upper edge
		const double* thickness {};		///< Cell's thickness along the
										///< normal (m)
		const double* n_R {};			///< Normal, R (or X) part
		const double* n_Z {};			///< Normal, Z (or Y) part
		const int* cartesian {};		///< 1 on a Cartesian zone
		const Params* params {};		///< nframes * n_elem
		std::size_t n_cells {};			///< Plasma cells next to the wall
		const std::uint64_t* cell_keys {};	///< Sorted cell_key values
		const int* cell_start {};		///< n_cells + 1 offsets
		const int* cell_elems {};		///< Elements of each cell
	};

	/**
	* Packs a cell's zone and indices into one key (each is far less than
	* 2^16).
	*/
	FLAN_HD inline std::uint64_t cell_key(const int zone, const int i,
		const int j, const int k)
	{
		auto u = [](const int v)
			{return static_cast<std::uint64_t>(v) & 0xffff;};
		return (u(zone) << 48) | (u(i) << 32) | (u(j) << 16) | u(k);
	}

	/**
	* Finds the wall element nearest a particle in a plasma cell next to
	* the wall, and the particle's distance from it: the fraction of the
	* way across the cell from the element's face, times the cell's
	* thickness along the element's normal. The first nearest in element
	* order is kept.
	*
	* @param s The sheath.
	* @param zone Zone of the particle's cell.
	* @param i x index of the cell.
	* @param j y index of the cell.
	* @param k z index of the cell.
	* @param pos The particle's computational position, 3 values.
	* @param[out] element The nearest element.
	* @param[out] d The distance from it along its normal (m), >= 0.
	* @return False if the cell has no wall element (nothing is set).
	*/
	FLAN_HD inline bool distance(const SheathView& s, const int zone,
		const int i, const int j, const int k, const double* pos,
		std::size_t& element, double& d)
	{
		if (s.n_elem == 0 || s.n_cells == 0) return false;
		const std::uint64_t key {cell_key(zone, i, j, k)};
		std::size_t left {0};
		std::size_t right {s.n_cells};
		while (left < right)
		{
			const std::size_t mid {(left + right) >> 1};
			if (s.cell_keys[mid] < key) left = mid + 1;
			else right = mid;
		}
		if (left >= s.n_cells || s.cell_keys[left] != key) return false;

		bool found {};
		for (int m {s.cell_start[left]}; m < s.cell_start[left + 1]; ++m)
		{
			const int e {s.cell_elems[m]};
			const int axis {s.face[e] / 2};
			const bool max_face {(s.face[e] & 1) == 1};
			const double frac {(max_face ? s.hi[e] - pos[axis]
				: pos[axis] - s.lo[e]) / (s.hi[e] - s.lo[e])};
			const double de {std::fmax(frac, 0.0) * s.thickness[e]};
			if (!found || de < d)
			{
				d = de;
				element = static_cast<std::size_t>(e);
				found = true;
			}
		}
		return found;
	}

	/// An element's Params in a frame, the frame held to the run's
	FLAN_HD inline const Params& params_at(const SheathView& s,
		const int frame, const std::size_t e)
	{
		const int t {frame < 0 ? 0 : (frame > s.nframes - 1 ? s.nframes - 1
			: frame)};
		return s.params[static_cast<std::size_t>(t) * s.n_elem + e];
	}

	/**
	* The sheath field at a particle (V/m, Cartesian), to add to the
	* background's: 0 unless the particle's cell is next to the wall, and
	* then the field of the nearest element along its normal (at the
	* particle's toroidal angle on a non-Cartesian zone).
	*
	* @param s The sheath.
	* @param zone Zone of the particle's cell.
	* @param frame Time frame (the particle's tidx).
	* @param i x index of the cell.
	* @param j y index of the cell.
	* @param k z index of the cell.
	* @param pos The particle's computational position, 3 values.
	* @param[out] E The field (X, Y, Z).
	*/
	FLAN_HD inline void field_at(const SheathView& s, const int zone,
		const int frame, const int i, const int j, const int k,
		const double* pos, double* E)
	{
		E[0] = 0.0;
		E[1] = 0.0;
		E[2] = 0.0;
		std::size_t e {};
		double d {};
		if (s.model == Model::off || !distance(s, zone, i, j, k, pos, e, d))
			return;
		const Params& p {params_at(s, frame, e)};
		if (d > range(p, s.model)) return;

		const double en {field(p, s.model, d)};
		if (s.cartesian[e])
		{
			E[0] = en * s.n_R[e];
			E[1] = en * s.n_Z[e];
			return;
		}
		const double phi {pos[2]};
		E[0] = en * s.n_R[e] * std::cos(phi);
		E[1] = en * s.n_R[e] * s.phi_sign * std::sin(phi);
		E[2] = en * s.n_Z[e];
	}

	/**
	* How many sub-steps a particle needs this time step: substeps for its
	* nearest wall element, or 1 away from the wall, with no charge, or
	* with the sheath off.
	*
	* @param s The sheath.
	* @param zone Zone of the particle's cell.
	* @param frame Time frame (the particle's tidx).
	* @param i x index of the cell.
	* @param j y index of the cell.
	* @param k z index of the cell.
	* @param pos The particle's computational position, 3 values.
	* @param speed Its speed (m/s).
	* @param q_over_m Its charge over mass (C/kg), > 0 for positive ions.
	* @param dt The time step (s).
	* @return The number of sub-steps, 1 if it needs none.
	*/
	FLAN_HD inline int substeps(const SheathView& s, const int zone,
		const int frame, const int i, const int j, const int k,
		const double* pos, const double speed, const double q_over_m,
		const double dt)
	{
		if (s.model == Model::off || !(q_over_m > 0.0)) return 1;
		std::size_t e {};
		double d {};
		if (!distance(s, zone, i, j, k, pos, e, d)) return 1;
		return substeps(params_at(s, frame, e), s.model, d, speed, q_over_m,
			dt);
	}

	/**
	* A sheath copied to a GPU (Sheath::Sheath::to_device): the view's
	* arrays are on that device. Freed with free_sheath.
	*/
	struct SheathDevice
	{
		int device_id {-1};		///< GPU the arrays are on, -1 for none
		SheathView view {};		///< The sheath, with device pointers
	};

}  // namespace Sheath
