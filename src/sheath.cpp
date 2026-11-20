/**
* @file sheath.cpp
*
* @brief The electric field of the sheath in front of a wall
*
* Between a plasma and a wall at an oblique magnetic field sit two thin
* layers: the magnetic presheath (Chodura layer), a few ion sound
* gyroradii thick, where ions are turned from along B toward the wall,
* and the Debye sheath, a few Debye lengths thick. Together they hold the
* floating sheath drop. They are far thinner than a background cell, so
* the background's potential doesn't have them; this adds their field as
* an analytic function of the distance from the wall, d, along the wall's
* normal. The collisional presheath (the 0.5 Te before the magnetic
* presheath) is in the background's own potential, so it isn't added.
*
* The drop is the floating one for no secondary emission,
*
*     lam_w = e phi_wall / Te = 0.5 ln[2 pi (m_e / m_i) (1 + Ti / Te)],
*
* about -2.84 for D at Ti = Te. The magnetic presheath takes
* lam_mps = ln cos(theta), theta the angle between B and the wall's normal
* (Stangeby 2012, Borodkina and Tsvetkov 2015), and the Debye sheath the
* rest, so the share in the Debye sheath is
*
*     f_D = 1 - ln cos(theta) / lam_w.
*
* f_D falls to 0 at Stangeby's critical angle, where the presheath takes
* the whole drop (about 86.6 degrees from the normal for D at Ti = Te),
* and is 0 past it: there is no Debye sheath there, and the flow reaches
* the wall subsonic. Kinetic simulations (Coulette and Manfredi 2016) find
* the same thing, but with a smooth change instead of a sharp angle. f_D
* itself goes to 0 smoothly (only its slope jumps), so it isn't held at a
* floor: the published results give no value for one.
*
* Two shapes for the potential (sheath_model):
*
* - "brooks": Brooks (1990), as in ERO and GITR,
*
*       phi(d) = phi_wall [f_D exp(-d / 2 lambda_D)
*                + (1 - f_D) exp(-d / rho_s)].
*
*   lambda_D is the electron Debye length sqrt(eps0 Te / (e n_e)) and
*   rho_s = c_s / Omega_i the ion sound gyroradius, with
*   c_s = sqrt(e (Te + Ti) / m_i). At Ti = Te that's GITR's
*   sqrt(2 kT / m) / Omega and Borodkina's rho_iCs; Brooks' own R_DT is
*   sqrt(kT / m) / Omega, sqrt(2) smaller.
*
* - "borodkina": Borodkina and Tsvetkov (2015). With xi = d / lambda_D,
*   a Debye layer from the wall to xi_mps,
*
*       lam(xi) = lam_w + Q (1 - exp(-a xi)),
*
*   then the presheath, lam(xi) = lam_mps exp(-(xi - xi_mps) / L), whose
*   decay length L lambda_D = rho_s sin(theta) grows with the angle. a, Q
*   and xi_mps are their eqs. (8)-(9), which make lam continuous at xi_mps
*   (and E too, to within their approximations). At normal incidence the
*   presheath is gone and the Debye layer is an exponential with
*   1 / a close to 2 lambda_D, as in Brooks. Past the critical angle there
*   is no Debye layer, and the presheath runs from the wall:
*   lam = lam_w exp(-xi / L). The field is continuous at xi_mps to well
*   within a percent at every angle below the critical one, though close
*   to it (within a degree or so) a and Q take odd values: the Debye layer
*   holds only a few percent of the drop there, so it matters little.
*
* Either way phi is negative and the field E = -dphi/dd n (n the unit
* normal into the plasma) points into the wall, pulling positive ions in.
* An ion of charge Z gains Z e |phi_wall| crossing both layers.
*
* The field is added only in the plasma cells next to the wall (see
* Sheath::Sheath), so a cell should be several rho_s thick for the field
* to have died away before its far side.
*/
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#include "constants.h"
#include "errors.h"
#include "sheath.h"

namespace Sheath
{
	/**
	* Works out the sheath in front of a wall from the plasma next to it.
	*
	* Te at or below 0 gives no sheath (phi_wall = 0). Without a field
	* (bmag = 0) there's no magnetic presheath: f_D = 1, as at normal
	* incidence.
	*
	* @param te Electron temperature (eV).
	* @param ti Ion temperature (eV).
	* @param ne Electron density (m^-3), > 0.
	* @param bmag Field strength (T).
	* @param cos_theta |b . n|, the cosine of the angle between B and the
	*        wall's normal.
	* @param elec_mass_amu Electron mass (amu).
	* @param ion_mass_amu Main ion mass (amu).
	* @return The sheath's Params.
	*/
	Params make_params(const double te, const double ti, const double ne,
		const double bmag, const double cos_theta, const double elec_mass_amu,
		const double ion_mass_amu)
	{
		Params p {};
		if (!(te > 0.0) || !(ne > 0.0)) return p;
		p.te = te;

		const double mi {ion_mass_amu * Constants::amu_to_kg};
		const double c_s {std::sqrt(Constants::ev_to_j * (te + std::max(ti,
			0.0)) / mi)};
		const double omega {Constants::ev_to_j * bmag / mi};
		p.lambda_d = std::sqrt(Constants::eps0 * te / (ne
			* Constants::ev_to_j));
		p.rho_s = omega > 0.0 ? c_s / omega
			: std::numeric_limits<double>::infinity();

		// Without a field there's no presheath, as at normal incidence
		const double c {omega > 0.0 ? std::clamp(std::abs(cos_theta), 0.0,
			1.0) : 1.0};
		p.theta = std::acos(c);
		const double s {std::sin(p.theta)};

		// The floating drop and its split
		p.lam_w = 0.5 * std::log(2.0 * Constants::pi * elec_mass_amu
			/ ion_mass_amu * (1.0 + std::max(ti, 0.0) / te));
		p.phi_wall = p.lam_w * te;
		const double ln_cos {c > 0.0 ? std::log(c)
			: -std::numeric_limits<double>::infinity()};
		p.lam_mps = std::max(ln_cos, p.lam_w);
		p.f_d = std::clamp(1.0 - p.lam_mps / p.lam_w, 0.0, 1.0);

		// Borodkina and Tsvetkov's profile. L is the presheath's decay
		// length in Debye lengths.
		p.l_mps = omega > 0.0 ? p.rho_s * s / p.lambda_d : 0.0;
		p.debye = ln_cos > p.lam_w;
		if (p.debye)
		{
			// -Delta n_mps lam_mps = (ln cos / sin)^2 (lambda_D / rho_s)^2,
			// written so that it goes to 0 at normal incidence
			const double g {s > 0.0 ? ln_cos / s : 0.0};
			const double ratio {omega > 0.0 ? p.lambda_d / p.rho_s : 0.0};
			const double dn_lam {g * g * ratio * ratio};
			const double drop {p.lam_w - p.lam_mps};
			const double c1 {dn_lam - 6.0 * c};
			const double root {2.0 * std::exp(p.lam_w) + 4.0 * c
				* std::sqrt(1.0 - drop) + c1};
			const double sq {std::sqrt(std::max(root, 0.0))};
			p.a = (std::sqrt(dn_lam) - sq) / drop;
			p.q = p.a != 0.0 ? sq / p.a : 0.0;
			// At normal incidence the presheath has no drop, Q = -lam_w and
			// the Debye layer goes on forever (xi_mps infinite)
			const double arg {(drop + p.q) / p.q};
			p.xi_mps = arg > 0.0 && ln_cos < 0.0 ? -std::log(arg) / p.a
				: std::numeric_limits<double>::infinity();

			// If the approximation stops working (a or Q not positive),
			// leave the Debye layer out
			if (!(p.a > 0.0) || !(p.q > 0.0) || !(p.xi_mps >= 0.0))
				p.debye = false;
		}
		if (!p.debye)
		{
			p.a = 0.0;
			p.q = 0.0;
			p.xi_mps = 0.0;
		}
		return p;
	}

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
	double potential(const Params& p, const Model model, const double d)
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
	double field(const Params& p, const Model model, const double d)
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
	double range(const Params& p, const Model model, const double n_scales)
	{
		if (model == Model::off || p.te == 0.0) return 0.0;
		const double l_ds {2.0 * p.lambda_d};
		if (model == Model::brooks)
		{
			return n_scales * (p.f_d < 1.0 ? std::max(l_ds, p.rho_s)
				: l_ds);
		}
		// Borodkina: the Debye layer then the presheath, or one of them
		if (p.debye && !std::isfinite(p.xi_mps))
			return n_scales / p.a * p.lambda_d;
		if (p.debye) return (p.xi_mps + n_scales * p.l_mps) * p.lambda_d;
		return n_scales * p.l_mps * p.lambda_d;
	}

	namespace
	{
		/// Packs a cell's zone and indices into one key (each is far less
		/// than 2^16)
		std::uint64_t cell_key(const int zone, const int i, const int j,
			const int k)
		{
			auto u = [](int v)
				{return static_cast<std::uint64_t>(v) & 0xffff;};
			return (u(zone) << 48) | (u(i) << 32) | (u(j) << 16) | u(k);
		}
	}

	/**
	* Works out the sheath in front of every wall element in every time
	* frame, and which elements each plasma cell next to the wall has.
	*
	* Each element's Params come from its plasma cell's Te, Ti, n_e and B,
	* with theta from B and the element's normal (in the middle of its phi
	* range), as Sputtering::Source does for the ion flux.
	*
	* A particle's distance from an element is measured across its cell:
	* the fraction of the way from the element's face to the opposite face
	* (in the computational coordinate across the face), times the cell's
	* thickness along the element's normal, taken as 2 (c - m) . n with c
	* the cell's center and m the middle of the face, in the poloidal
	* plane (the X-Y plane on a Cartesian zone). On a Cartesian zone, and
	* on any cell that is a rectangle in that plane, that's the exact
	* distance. Elements whose cell is thinner than 7 of the field's
	* longest scale lengths in some frame are counted (num_thin): the field
	* is cut off at the cell's far side there while still about 1e-3 of
	* its size at the wall.
	*
	* @param surface The wall surface elements.
	* @param bkg The background, with its plasma, field and geometry.
	* @param model The profile. With Model::off nothing is set up.
	*/
	Sheath::Sheath(const Sputtering::Surface& surface,
		const Background::Background& bkg, const Model model)
		: m_model {model},
		m_nframes {static_cast<int>(bkg.get_times().size())},
		m_phi_sign {bkg.get_phi_sign()}
	{
		if (model == Model::off) return;

		const Zones::Zones& zones {bkg.get_zones()};
		const auto& X {bkg.get_X()};
		const auto& Y {bkg.get_Y()};
		const auto& Z {bkg.get_Z()};
		const auto& ne {bkg.get_ne()};
		const auto& te {bkg.get_te()};
		const auto& ti {bkg.get_ti()};
		const auto& bX {bkg.get_bX()};
		const auto& bY {bkg.get_bY()};
		const auto& bZ {bkg.get_bZ()};
		const double me {bkg.get_elec_mass_amu()};
		const double mi {bkg.get_ion_mass_amu()};

		const auto& els {surface.elements()};
		const std::size_t n_elem {els.size()};
		m_elements.reserve(n_elem);
		std::vector<int> thin(n_elem, 0);
		for (std::size_t e {}; e < n_elem; ++e)
		{
			const auto& el {els[e]};
			Element s {};
			s.zone = el.zone;
			s.face = el.face;
			s.n_R = el.n_R;
			s.n_Z = el.n_Z;
			s.cartesian = el.cartesian;

			// The cell's extent across the face
			const int axis {static_cast<int>(el.face) / 2};
			const auto& grid {axis == 0 ? zones.get_grid_x(el.zone)
				: zones.get_grid_y(el.zone)};
			const int idx {axis == 0 ? el.i : el.j};
			s.lo = grid[idx];
			s.hi = grid[idx + 1];

			// Its center in the plane of the element's (R, Z) or (X, Y)
			double cR {};
			double cZ {};
			if (el.cartesian)
			{
				const auto& gx {zones.get_grid_x(el.zone)};
				const auto& gy {zones.get_grid_y(el.zone)};
				cR = 0.5 * (gx[el.i] + gx[el.i + 1]);
				cZ = 0.5 * (gy[el.j] + gy[el.j + 1]);
			}
			else
			{
				cR = std::hypot(X(el.zone, el.i, el.j, el.k),
					Y(el.zone, el.i, el.j, el.k));
				cZ = Z(el.zone, el.i, el.j, el.k);
			}
			s.thickness = 2.0 * std::abs((cR - el.R_mid) * el.n_R
				+ (cZ - el.Z_mid) * el.n_Z);
			if (!(s.thickness > 0.0))
			{
				Errors::fatal_error("Sheath: the plasma cell of a wall "
					"element has no thickness along the element's "
					"normal.");
			}
			m_elements.push_back(s);
			m_cells[cell_key(el.zone, el.i, el.j, el.k)].push_back(e);
		}

		// Params per frame and element
		m_params.resize(static_cast<std::size_t>(m_nframes) * n_elem);
		for (int t {}; t < m_nframes; ++t)
		for (std::size_t e {}; e < n_elem; ++e)
		{
			const auto& el {els[e]};
			const int z {el.zone};
			const auto n {el.normal_xyz(0.5 * (el.phi_min + el.phi_max),
				m_phi_sign)};
			const std::array<double, 3> B {bX(z, t, el.i, el.j, el.k),
				bY(z, t, el.i, el.j, el.k), bZ(z, t, el.i, el.j, el.k)};
			const double bmag {std::sqrt(B[0] * B[0] + B[1] * B[1]
				+ B[2] * B[2])};
			const double cos_theta {bmag > 0.0
				? (B[0] * n[0] + B[1] * n[1] + B[2] * n[2]) / bmag : 1.0};
			Params& p {m_params[t * n_elem + e]};
			p = make_params(te(z, t, el.i, el.j, el.k),
				ti(z, t, el.i, el.j, el.k), ne(z, t, el.i, el.j, el.k), bmag,
				cos_theta, me, mi);
			// Thin: the field is still above 1e-3 of the drop per scale
			// length at the cell's far side
			if (range(p, m_model, 7.0) > m_elements[e].thickness)
			{
				thin[e] = 1;
				if (p.phi_wall < 0.0)
				{
					m_max_lost = std::max(m_max_lost, potential(p, m_model,
						m_elements[e].thickness) / p.phi_wall);
				}
			}
		}
		for (const int v : thin) m_num_thin += v;
	}

	/**
	* @name Accessors
	* The profile, the number of elements and frames, the Params of an
	* element in a frame, how many elements have a cell thinner than the
	* field's reach, and the largest fraction of the drop beyond a cell's
	* far side, which those particles don't gain.
	*/
	///@{
	Model Sheath::model() const {return m_model;}
	std::size_t Sheath::size() const {return m_elements.size();}
	int Sheath::nframes() const {return m_nframes;}
	const Params& Sheath::params(const int frame,
		const std::size_t element) const
	{
		return m_params[static_cast<std::size_t>(frame) * m_elements.size()
			+ element];
	}
	int Sheath::num_thin() const {return m_num_thin;}
	double Sheath::max_lost() const {return m_max_lost;}
	///@}

	/**
	* Finds the wall element nearest a particle in a plasma cell next to
	* the wall, and the particle's distance from it.
	*
	* @param zone Zone of the particle's cell.
	* @param i x index of the cell.
	* @param j y index of the cell.
	* @param k z index of the cell.
	* @param pos The particle's computational position (x, y, z).
	* @param element Set to the nearest element.
	* @param d Set to the distance from it along its normal (m), >= 0.
	* @return False if the cell has no wall element (nothing is set).
	*/
	bool Sheath::distance(const int zone, const int i, const int j,
		const int k, const std::array<double, 3>& pos,
		std::size_t& element, double& d) const
	{
		if (m_elements.empty()) return false;
		const auto it {m_cells.find(cell_key(zone, i, j, k))};
		if (it == m_cells.end()) return false;

		bool found {};
		for (const std::size_t e : it->second)
		{
			const Element& s {m_elements[e]};
			const int axis {static_cast<int>(s.face) / 2};
			const bool max_face {(static_cast<int>(s.face) & 1) == 1};
			const double frac {(max_face ? s.hi - pos[axis]
				: pos[axis] - s.lo) / (s.hi - s.lo)};
			const double de {std::max(frac, 0.0) * s.thickness};
			if (!found || de < d)
			{
				d = de;
				element = e;
				found = true;
			}
		}
		return found;
	}

	/**
	* The sheath field at a particle (V/m, Cartesian), to add to the
	* background's. It is 0 unless the particle's cell is next to the
	* wall; then it's the field of the nearest element along its normal.
	*
	* @param zone Zone of the particle's cell.
	* @param frame Time frame (the particle's tidx).
	* @param i x index of the cell.
	* @param j y index of the cell.
	* @param k z index of the cell.
	* @param pos The particle's computational position (x, y, z).
	* @return E (X, Y, Z).
	*/
	std::array<double, 3> Sheath::field_at(const int zone, const int frame,
		const int i, const int j, const int k,
		const std::array<double, 3>& pos) const
	{
		std::size_t e {};
		double d {};
		if (m_model == Model::off || !distance(zone, i, j, k, pos, e, d))
			return {0.0, 0.0, 0.0};
		const int t {std::clamp(frame, 0, m_nframes - 1)};
		const Params& p {params(t, e)};
		if (d > range(p, m_model)) return {0.0, 0.0, 0.0};

		const Element& s {m_elements[e]};
		const double en {field(p, m_model, d)};
		if (s.cartesian) return {en * s.n_R, en * s.n_Z, 0.0};
		const double phi {pos[2]};
		return {en * s.n_R * std::cos(phi),
			en * s.n_R * m_phi_sign * std::sin(phi), en * s.n_Z};
	}

	/**
	* How many sub-steps a particle needs this time step so that none moves
	* it more than substep_frac of the field's local scale length.
	*
	* A particle that can't reach the field this step (or has no charge)
	* needs 1, a normal step; one that can, at least 2. Otherwise the scale length is the Debye
	* layer's (2 lambda_D for "brooks", lambda_D / a for "borodkina") if the
	* particle can reach the Debye layer this step, and the presheath's
	* (rho_s, or rho_s sin(theta)) if not. The speed it might reach is its
	* own plus what the rest of the drop could add,
	* sqrt(2 (q/m) (phi(d) - phi_wall)), so the count holds as it speeds
	* up toward the wall. The count is capped at max_substeps.
	*
	* @param zone Zone of the particle's cell.
	* @param frame Time frame (the particle's tidx).
	* @param i x index of the cell.
	* @param j y index of the cell.
	* @param k z index of the cell.
	* @param pos The particle's computational position (x, y, z).
	* @param speed The particle's speed (m/s).
	* @param q_over_m Its charge over mass (C/kg), > 0 for positive ions.
	* @param dt The time step (s).
	* @return The number of sub-steps, 1 if it needs none.
	*/
	int Sheath::substeps(const int zone, const int frame, const int i,
		const int j, const int k, const std::array<double, 3>& pos,
		const double speed, const double q_over_m, const double dt) const
	{
		if (m_model == Model::off || !(q_over_m > 0.0)) return 1;
		std::size_t e {};
		double d {};
		if (!distance(zone, i, j, k, pos, e, d)) return 1;
		const int t {std::clamp(frame, 0, m_nframes - 1)};
		const Params& p {params(t, e)};
		if (p.te == 0.0) return 1;

		// The fastest it could be going, and how close it could get
		const double drop {std::max(potential(p, m_model, d) - p.phi_wall,
			0.0)};
		const double v_max {speed + std::sqrt(2.0 * q_over_m * drop)};
		const double reach {v_max * dt};
		const double nearest {d - reach};
		if (nearest > range(p, m_model, 7.0)) return 1;

		// The Debye layer's scale and extent, and the presheath's scale
		double l_debye {};
		double debye_edge {};
		double l_mps {};
		if (m_model == Model::brooks)
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
			scale = std::min(l_debye, l_mps);
		if (!(scale > 0.0)) return 1;

		// At least 2 in the field, so that a particle crossing the wall
		// gets its impact velocity (see Boris::substep_cpu)
		const double n {std::ceil(reach / (substep_frac * scale))};
		return static_cast<int>(std::clamp(n, 2.0,
			static_cast<double>(max_substeps)));
	}

	/**
	* Prints the sheath model, the range of the drop, of f_D and of the
	* field angle over the elements and frames, and how many elements have
	* a cell thinner than the field's reach.
	*
	* @param os Stream to print to.
	*/
	void Sheath::print_summary(std::ostream& os) const
	{
		if (m_model == Model::off) return;
		os << "Sheath field (" << (m_model == Model::brooks ? "brooks"
			: "borodkina") << ") on " << m_elements.size() << " elements\n";
		if (m_params.empty()) return;
		double drop_lo {1e300}, drop_hi {}, fd_lo {1.0}, fd_hi {};
		double th_lo {1e300}, th_hi {};
		for (const Params& p : m_params)
		{
			if (p.te == 0.0) continue;
			drop_lo = std::min(drop_lo, -p.phi_wall);
			drop_hi = std::max(drop_hi, -p.phi_wall);
			fd_lo = std::min(fd_lo, p.f_d);
			fd_hi = std::max(fd_hi, p.f_d);
			th_lo = std::min(th_lo, p.theta);
			th_hi = std::max(th_hi, p.theta);
		}
		const double deg {180.0 / Constants::pi};
		os << "  Drop " << drop_lo << " to " << drop_hi << " V, f_D "
			<< fd_lo << " to " << fd_hi << ", field " << th_lo * deg
			<< " to " << th_hi * deg << " degrees from the normal\n";
		if (m_num_thin > 0)
		{
			os << "  Warning: " << m_num_thin << " elements have a plasma "
				"cell thinner than the sheath's reach, so the field is cut "
				"off at the cell's far side, leaving out up to "
				<< 100.0 * m_max_lost << "% of the drop\n";
		}
	}
}
