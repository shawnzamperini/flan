#pragma once
/**
* @file sputtering_hd.h
*
* @brief Starting sputtered atoms and tallying where they land, shared by
*        the CPU and the GPU (see host_device.h).
*
* Sputtering::Source::sample_birth and Sputtering::Deposition's
* add_into_wall_cell and add_at_wall_face (sputtering.cpp) call these on
* views of their own data; the GPU calls them on device copies of the same
* data (SourceDevice, DepositionDevice). So a GPU run starts atoms and
* tallies their landings the same way as a CPU run.
*
* Views hold plain pointers and sizes, never std::vector, so they can be
* copied to a kernel as they are. The tallies are added to with an adder
* given as a template parameter: the CPU's adds with omp atomic, the GPU's
* with atomicAdd. (An inline function containing an omp pragma would be
* compiled differently by nvcc and the host compiler, and the linker could
* keep either.)
*/

#include <cmath>
#include <cstddef>
#include <cstdint>

#include "constants.h"
#include "host_device.h"
#include "pcg32.h"
#include "zones_hd.h"


namespace Sputtering
{
	/**
	* The surface elements as flat arrays, one value per element (see
	* SurfaceElement for what each is).
	*/
	struct ElementsView
	{
		std::size_t n {};				///< Number of elements
		const int* zone {};				///< Zone of the plasma cell
		const int* i {};				///< x index of the plasma cell
		const int* j {};				///< y index of the plasma cell
		const int* k {};				///< z index of the plasma cell
		const int* face {};				///< Zones::Face on the wall
		const double* n_R {};			///< R part of the normal
		const double* n_Z {};			///< Z part of the normal
		const double* phi_min {};		///< Toroidal angle at one side
		const double* phi_max {};		///< At the other side
		const int* cartesian {};		///< 1 on a Cartesian zone
	};

	/**
	* Unit normal of an element into the plasma, in Cartesian (X, Y, Z),
	* at toroidal angle z (SurfaceElement::normal_xyz).
	*
	* @param el The elements.
	* @param e Index of the element.
	* @param z Toroidal angle (rad), unused on a Cartesian zone.
	* @param phi_sign Handedness of (R, phi, Z), +1 or -1.
	* @param[out] n The normal, 3 values.
	*/
	FLAN_HD inline void element_normal(const ElementsView& el,
		const std::size_t e, const double z, const double phi_sign, double* n)
	{
		if (el.cartesian[e])
		{
			n[0] = el.n_R[e];
			n[1] = el.n_Z[e];
			n[2] = 0.0;
			return;
		}
		n[0] = el.n_R[e] * std::cos(z);
		n[1] = el.n_R[e] * phi_sign * std::sin(z);
		n[2] = el.n_Z[e];
	}

	/**
	* Lower index of the interval of a grid that holds a value: one less
	* than the first grid value greater than it (std::upper_bound), held to
	* 0 to n - 2.
	*
	* @param grid Strictly increasing values, at least two.
	* @param n Number of values.
	* @param value Value to find.
	* @return Index i, with value between grid[i] and grid[i + 1] if it's
	*         inside the grid.
	*/
	FLAN_HD inline int lower_index(const double* grid, const int n,
		const double value)
	{
		int left {0};
		int right {n};
		while (left < right)
		{
			const int mid {(left + right) >> 1};
			if (grid[mid] <= value) left = mid + 1;
			else right = mid;
		}
		const int i {left - 1 > 0 ? left - 1 : 0};
		return i < n - 2 ? i : n - 2;
	}

	/**
	* A sputtering yield table as flat arrays (see YieldTable).
	*/
	struct YieldView
	{
		int n_energy {};				///< Number of energies
		int n_angle {};					///< Number of angles
		const double* log_energy {};	///< log10(energy [eV])
		const double* angle {};			///< Angle from the normal (deg)
		const double* yield {};			///< (energy, angle), row major
		double energy_min {};			///< Lowest energy (eV)
		double target_Es {};			///< Surface binding energy (eV)
	};

	/**
	* Yield at an energy and angle (YieldTable::get_yield): bilinear in
	* log10(energy) and angle, 0 below the table's lowest energy, the
	* highest energy's yield above it, and the angle held to the table.
	*
	* @param y The table.
	* @param energy Energy of the projectile (eV).
	* @param angle Angle from the normal (degrees).
	* @return The yield.
	*/
	FLAN_HD inline double get_yield(const YieldView& y, const double energy,
		double angle)
	{
		if (!(energy >= y.energy_min)) return 0.0;

		const double log_energy {std::fmin(std::log10(energy),
			y.log_energy[y.n_energy - 1])};
		angle = std::fmin(std::fmax(angle, y.angle[0]),
			y.angle[y.n_angle - 1]);

		const int i {lower_index(y.log_energy, y.n_energy, log_energy)};
		const int j {lower_index(y.angle, y.n_angle, angle)};
		const double fe {(log_energy - y.log_energy[i])
			/ (y.log_energy[i + 1] - y.log_energy[i])};
		const double fa {(angle - y.angle[j])
			/ (y.angle[j + 1] - y.angle[j])};

		const int na {y.n_angle};
		const double y00 {y.yield[i * na + j]};
		const double y01 {y.yield[i * na + j + 1]};
		const double y10 {y.yield[(i + 1) * na + j]};
		const double y11 {y.yield[(i + 1) * na + j + 1]};

		return (1.0 - fe) * ((1.0 - fa) * y00 + fa * y01)
			+ fe * ((1.0 - fa) * y10 + fa * y11);
	}

	/**
	* Nodes x and weights w for averaging over the thermal energies of
	* the ions hitting a wall (see thermal_quadrature in sputtering.cpp).
	*/
	struct Quadrature
	{
		int n {};				///< Number of nodes
		const double* x {};		///< Energy over Ti
		const double* w {};		///< Weight, summing to 1
	};

	/**
	* Samples the energy of the ion that sputtered an atom, for the
	* formula: e_sheath plus a thermal energy from Gamma(2, ti), weighted
	* by the yield at each energy (see sample_impact_energy in
	* sputtering.cpp). The quadrature node is picked by its share of
	* w Y, the first whose running sum is above u times the total.
	*
	* The running sum is formed twice, the same way, rather than stored,
	* so nothing is allocated (the GPU can't).
	*
	* @param y The yield table.
	* @param q The thermal quadrature.
	* @param e_sheath Energy gained in the sheath and presheath (eV).
	* @param ti Ion temperature (eV).
	* @param angle Angle of incidence (degrees).
	* @param u Uniform random number in [0, 1).
	* @return The ion's energy (eV). e_sheath + 2 ti (the mean) if
	*         nothing is sputtered at any energy.
	*/
	FLAN_HD inline double sample_impact_energy(const YieldView& y,
		const Quadrature& q, const double e_sheath, const double ti,
		const double angle, const double u)
	{
		double total {};
		for (int n {}; n < q.n; ++n)
			total += q.w[n] * get_yield(y, e_sheath + q.x[n] * ti, angle);
		if (!(total > 0.0)) return e_sheath + 2.0 * ti;

		const double target {u * total};
		double running {};
		int pick {q.n - 1};
		for (int n {}; n < q.n; ++n)
		{
			running += q.w[n] * get_yield(y, e_sheath + q.x[n] * ti, angle);
			if (running > target)
			{
				pick = n;
				break;
			}
		}
		return e_sheath + q.x[pick] * ti;
	}

	/**
	* Picks the energy of a sputtered atom from the Thompson distribution,
	* cut off at e_max.
	*
	* The Thompson distribution is f(E) ~ E / (E + e_b)^3, with e_b the
	* surface binding energy. Its cumulative distribution is
	* F(E) = E^2 / (E + e_b)^2, so with u' = u F(e_max), solving F(E) = u'
	* gives E = e_b sqrt(u') / (1 - sqrt(u')), which is exactly distributed
	* as f(E) between 0 and e_max.
	*
	* @param e_b Surface binding energy (eV), > 0.
	* @param e_max Largest energy (eV). If it isn't > 0, e_b is returned.
	* @param u Uniform random number in [0, 1).
	* @return Energy of the sputtered atom (eV).
	*/
	FLAN_HD inline double sample_thompson(const double e_b,
		const double e_max, const double u)
	{
		if (!(e_max > 0.0)) return e_b;
		const double f_max {e_max * e_max / ((e_max + e_b) * (e_max + e_b))};
		const double s {std::sqrt(u * f_max)};
		return e_b * s / (1.0 - s);
	}

	/**
	* A direction from a cosine distribution about a unit normal (see
	* cosine_direction in sputtering.cpp).
	*
	* @param n The unit normal, 3 values.
	* @param u1 Uniform random number for the polar angle.
	* @param u2 Uniform random number for the azimuth.
	* @param[out] d The unit direction, 3 values.
	*/
	FLAN_HD inline void cosine_direction(const double* n, const double u1,
		const double u2, double* d)
	{
		// A perpendicular to n: cross it with the axis it's least along
		double axis[3] {0.0, 0.0, 0.0};
		const int least {std::fabs(n[0]) <= std::fabs(n[1])
			? (std::fabs(n[0]) <= std::fabs(n[2]) ? 0 : 2)
			: (std::fabs(n[1]) <= std::fabs(n[2]) ? 1 : 2)};
		axis[least] = 1.0;
		double t1[3] {n[1] * axis[2] - n[2] * axis[1],
			n[2] * axis[0] - n[0] * axis[2], n[0] * axis[1] - n[1] * axis[0]};
		const double len {std::sqrt(t1[0] * t1[0] + t1[1] * t1[1]
			+ t1[2] * t1[2])};
		for (int c {}; c < 3; ++c) t1[c] /= len;
		const double t2[3] {n[1] * t1[2] - n[2] * t1[1],
			n[2] * t1[0] - n[0] * t1[2], n[0] * t1[1] - n[1] * t1[0]};

		const double cos_th {std::sqrt(u1)};
		const double sin_th {std::sqrt(1.0 - u1)};
		const double az {2.0 * Constants::pi * u2};
		for (int c {}; c < 3; ++c)
		{
			d[c] = cos_th * n[c] + sin_th * (std::cos(az) * t1[c]
				+ std::sin(az) * t2[c]);
		}
	}

	/**
	* Where, when and how a sputtered atom starts (Source::sample_birth).
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
		double a_impact {};		///< Its angle from the normal (degrees)
		double e_atom {};		///< Energy of the sputtered atom (eV)
	};

	/**
	* What sample_birth needs from a Source, as flat arrays. idx is
	* frame * (number of elements) + element.
	*/
	struct SourceView
	{
		ElementsView el {};				///< The surface elements
		YieldView table {};				///< The yield table
		Quadrature quad {};				///< Thermal quadrature
		std::size_t n_idx {};			///< nframes * elements
		const double* cdf {};			///< Cumulative atoms/s per idx
		const double* frame_lo {};		///< Start of each frame's window
		const double* frame_hi {};		///< End of each frame's window
		const double* cell_bounds {};	///< x_lo, x_hi, y_lo, y_hi per
										///< element
		const double* e_sheath {};		///< Sheath energy (eV) per idx
		const double* ti {};			///< Ion temperature (eV) per idx
		const double* te {};			///< Electron temperature per idx
		const double* energy {};		///< Mean impact energy per idx
		const double* impact_angle {};	///< Mean impact angle per idx
		const int* point {};			///< Impact table point per idx, -1
										///< for the formula
		const int* sample_offset {};	///< First sample of each point
										///< (points + 1 values)
		const double* sample_energy {};	///< Each sample's energy / Te
		const double* sample_angle {};	///< Its angle (degrees)
		const double* sample_weight {};	///< Its share of the flux
		double angle {};				///< sput_angle_val (degrees)
		double lambda {};				///< Largest energy fraction an ion
										///< can give a target atom
		double imp_mass_kg {};			///< Mass of the sputtered atoms
		double phi_sign {1.0};			///< Handedness of (R, phi, Z)
	};

	/**
	* Starts a sputtered atom (see Source::sample_birth for the steps):
	* the frame and element in proportion to the atoms they sputter, the
	* time in the frame's window, the place along the face just inside the
	* plasma cell, the energy and angle of the ion that sputtered it, the
	* atom's Thompson energy and a cosine direction about the normal. The
	* random numbers are drawn in the same order on the CPU and the GPU.
	*
	* The source must sputter something (cdf[n_idx - 1] > 0); the CPU
	* checks this before calling.
	*
	* @param s The source.
	* @param rng Random number generator for this particle or thread.
	* @return The birth.
	*/
	FLAN_HD inline Birth sample_birth(const SourceView& s, pcg32& rng)
	{
		// 1. Frame and element, in proportion to the atoms they sputter:
		// the first idx whose cumulative rate is above r (upper_bound)
		const double r {rng.next_double() * s.cdf[s.n_idx - 1]};
		std::size_t left {0};
		std::size_t right {s.n_idx};
		while (left < right)
		{
			const std::size_t mid {(left + right) >> 1};
			if (s.cdf[mid] <= r) left = mid + 1;
			else right = mid;
		}
		const std::size_t idx {left < s.n_idx - 1 ? left : s.n_idx - 1};
		const std::size_t n_elem {s.el.n};
		Birth b {};
		b.frame = static_cast<int>(idx / n_elem);
		b.element = idx % n_elem;
		const std::size_t e {b.element};
		b.zone = s.el.zone[e];

		// Time, uniform in the frame's window
		b.t = s.frame_lo[b.frame] + rng.next_double()
			* (s.frame_hi[b.frame] - s.frame_lo[b.frame]);

		// 2. Place: uniform along the face, just inside the plasma cell
		constexpr double inset {1e-3};
		const double x_lo {s.cell_bounds[4 * e]};
		const double x_hi {s.cell_bounds[4 * e + 1]};
		const double y_lo {s.cell_bounds[4 * e + 2]};
		const double y_hi {s.cell_bounds[4 * e + 3]};
		const double ux {rng.next_double()};
		switch (s.el.face[e])
		{
			case Zones::xmin:
				b.x = x_lo + inset * (x_hi - x_lo);
				b.y = y_lo + ux * (y_hi - y_lo);
				break;
			case Zones::xmax:
				b.x = x_hi - inset * (x_hi - x_lo);
				b.y = y_lo + ux * (y_hi - y_lo);
				break;
			case Zones::ymin:
				b.x = x_lo + ux * (x_hi - x_lo);
				b.y = y_lo + inset * (y_hi - y_lo);
				break;
			default:	// ymax
				b.x = x_lo + ux * (x_hi - x_lo);
				b.y = y_hi - inset * (y_hi - y_lo);
				break;
		}
		b.z = s.el.phi_min[e] + rng.next_double()
			* (s.el.phi_max[e] - s.el.phi_min[e]);

		// 3. Energy (and angle) of the ion that sputtered it: from the
		// impact samples of the element's point, weighted by their share
		// of the flux times their yield, or from the formula at
		// sput_angle_val. The running sum over the samples is formed
		// twice, as in sample_impact_energy.
		const int point {s.point[idx]};
		if (point >= 0)
		{
			const int k0 {s.sample_offset[point]};
			const int k1 {s.sample_offset[point + 1]};
			auto share = [&](const int k)
			{
				return s.sample_weight[k] * get_yield(s.table,
					s.sample_energy[k] * s.te[idx], s.sample_angle[k]);
			};
			double total {};
			for (int k {k0}; k < k1; ++k) total += share(k);
			if (total > 0.0)
			{
				const double target {rng.next_double() * total};
				double running {};
				int pick {k1 - 1};
				for (int k {k0}; k < k1; ++k)
				{
					running += share(k);
					if (running > target)
					{
						pick = k;
						break;
					}
				}
				b.e_impact = s.sample_energy[pick] * s.te[idx];
				b.a_impact = s.sample_angle[pick];
			}
			else
			{
				b.e_impact = s.energy[idx];
				b.a_impact = s.impact_angle[idx];
			}
		}
		else
		{
			b.e_impact = sample_impact_energy(s.table, s.quad,
				s.e_sheath[idx], s.ti[idx], s.angle, rng.next_double());
			b.a_impact = s.angle;
		}

		// 4. Energy of the atom, Thompson cut off at Lambda E_ion
		b.e_atom = sample_thompson(s.table.target_Es,
			s.lambda * b.e_impact, rng.next_double());
		const double speed {std::sqrt(2.0 * b.e_atom * Constants::ev_to_j
			/ s.imp_mass_kg)};

		// 5. Direction: cosine about the normal at the atom's phi
		double n[3];
		element_normal(s.el, e, b.z, s.phi_sign, n);
		const double u1 {rng.next_double()};
		const double u2 {rng.next_double()};
		double d[3];
		cosine_direction(n, u1, u2, d);
		b.vX = speed * d[0];
		b.vY = speed * d[1];
		b.vZ = speed * d[2];
		return b;
	}

	/**
	* Packs a cell's zone and indices into one key (each is far less than
	* 2^16), for finding the elements on a wall cell.
	*/
	FLAN_HD inline std::uint64_t cell_key(const int zone, const int i,
		const int j, const int k)
	{
		auto u = [](const int v)
			{return static_cast<std::uint64_t>(v) & 0xffff;};
		return (u(zone) << 48) | (u(i) << 32) | (u(j) << 16) | u(k);
	}

	/**
	* Packs the cell of a zone's wall face an element is on into one key:
	* the zone, the face, the cell's index along the face (j for an x
	* face, i for a y face) and its toroidal cell. It has its own lookup,
	* so it can't clash with cell_key.
	*/
	FLAN_HD inline std::uint64_t wall_face_key(const int zone, const int face,
		const int i, const int j, const int k)
	{
		const bool x_face {face == Zones::xmin || face == Zones::xmax};
		return cell_key(zone, face, x_face ? j : i, k);
	}

	/**
	* Position of a key in a sorted array of keys, or -1 if it isn't
	* there.
	*/
	FLAN_HD inline long find_key(const std::uint64_t* keys,
		const std::size_t n, const std::uint64_t key)
	{
		std::size_t left {0};
		std::size_t right {n};
		while (left < right)
		{
			const std::size_t mid {(left + right) >> 1};
			if (keys[mid] < key) left = mid + 1;
			else right = mid;
		}
		return left < n && keys[left] == key ? static_cast<long>(left) : -1;
	}

	/// Indices of Deposition's scalar tallies in DepositionView::scalars
	enum DepositionScalar
	{
		dep_unassigned,		///< Into wall cells with no element
		dep_wall_face,		///< At wall faces with no element
		dep_lifetime_sum,	///< Sum of weight * time since birth
		dep_energy_sum,		///< Sum of weight * impact energy
		dep_energy_weight,	///< Rate with an impact energy
		dep_angle_sum,		///< Sum of weight * impact angle
		dep_angle_weight,	///< Rate with an impact angle
		num_dep_scalars
	};

	/**
	* Deposition's lookups and tallies as flat arrays. The tallies are
	* rates (atoms/s), added to as atoms land.
	*/
	struct DepositionView
	{
		ElementsView el {};				///< The surface elements
		double phi_sign {1.0};			///< Handedness of (R, phi, Z)

		/// Wall cells bordering the plasma, sorted by cell_key, and the
		/// elements on each (wall_cell_start[c] to wall_cell_start[c + 1]
		/// in wall_cell_elems, in element order)
		std::size_t n_wall_cells {};
		const std::uint64_t* wall_cell_keys {};
		const int* wall_cell_start {};
		const int* wall_cell_elems {};

		/// Cells of zone wall faces with an element, sorted by
		/// wall_face_key, and the element on each
		std::size_t n_wall_faces {};
		const std::uint64_t* wall_face_keys {};
		const int* wall_face_elems {};

		int n_lifetime {};				///< Time since birth bins
		const double* lifetime_edges {};	///< n_lifetime + 1 edges (s)
		int n_energy {};				///< Impact energy bins
		const double* energy_edges {};	///< n_energy + 1 edges (eV)
		int n_angle {};					///< 1 degree angle bins
		int n_charge {};				///< Charge states, 0 to Z

		double* redeposited {};		///< Rate onto each element
		double* prompt {};			///< Prompt rate onto each element
		double* element_energy {};	///< Sum of weight * impact energy on
									///< each element
		double* lifetime {};		///< Rate in each time bin
		double* charge {};			///< Rate at each charge state
		double* prompt_charge {};	///< Prompt rate at each charge state
		double* energy {};			///< Rate in each energy bin
		double* angle {};			///< Rate in each angle bin
		double* scalars {};			///< num_dep_scalars values
	};

	/**
	* Bin of a value in n bins with n + 1 increasing edges (upper_bound
	* minus 1), with values outside them in the first or last bin.
	*/
	FLAN_HD inline int edge_bin(const double* edges, const int n,
		const double value)
	{
		int left {0};
		int right {n + 1};
		while (left < right)
		{
			const int mid {(left + right) >> 1};
			if (edges[mid] <= value) left = mid + 1;
			else right = mid;
		}
		if (left == 0) return 0;
		return left - 1 < n - 1 ? left - 1 : n - 1;
	}

	/**
	* Adds an absorbed particle's impact energy and angle to the tallies
	* (see Deposition::add_impact): the energy always, the angle from the
	* element's normal at toroidal angle z if it landed on element e.
	*
	* @param d The tallies.
	* @param v Its velocity (vX, vY, vZ) (m/s), or nullptr to add nothing.
	* @param mass Its mass (kg).
	* @param e Element it landed on, or -1.
	* @param z Its toroidal angle, for the normal.
	* @param weight Its weight (atoms/s).
	* @param add Adder, add(double* p, double value).
	*/
	template <typename Add>
	FLAN_HD inline void add_impact(const DepositionView& d, const double* v,
		const double mass, const long e, const double z, const double weight,
		const Add& add)
	{
		if (v == nullptr) return;
		const double v2 {v[0] * v[0] + v[1] * v[1] + v[2] * v[2]};
		const double energy {0.5 * mass * v2 / Constants::ev_to_j};
		add(&d.energy[edge_bin(d.energy_edges, d.n_energy, energy)], weight);
		add(&d.scalars[dep_energy_sum], weight * energy);
		add(&d.scalars[dep_energy_weight], weight);

		if (e < 0 || !(v2 > 0.0)) return;
		add(&d.element_energy[e], weight * energy);

		// The normal points into the plasma, the velocity into the wall
		double n[3];
		element_normal(d.el, static_cast<std::size_t>(e), z, d.phi_sign, n);
		const double c {std::fmin(1.0, std::fmax(-1.0, -(v[0] * n[0]
			+ v[1] * n[1] + v[2] * n[2]) / std::sqrt(v2)))};
		const double angle {std::acos(std::fabs(c)) * 180.0
			/ Constants::pi};
		int abin {static_cast<int>(angle)};
		abin = abin < 0 ? 0 : (abin > d.n_angle - 1 ? d.n_angle - 1 : abin);
		add(&d.angle[abin], weight);
		add(&d.scalars[dep_angle_sum], weight * angle);
		add(&d.scalars[dep_angle_weight], weight);
	}

	/**
	* Adds an absorbed atom to the time since birth and charge tallies,
	* and to the prompt charge tally if it came back promptly (see
	* Deposition::add_lifetime_and_charge).
	*
	* @return True if the return was prompt.
	*/
	template <typename Add>
	FLAN_HD inline bool add_lifetime_and_charge(const DepositionView& d,
		const double pt, const int q, const double weight,
		const double prompt_by, const Add& add)
	{
		add(&d.lifetime[edge_bin(d.lifetime_edges, d.n_lifetime, pt)],
			weight);
		add(&d.scalars[dep_lifetime_sum], weight * pt);
		const int charge {q < 0 ? 0 : (q > d.n_charge - 1 ? d.n_charge - 1
			: q)};
		add(&d.charge[charge], weight);
		const bool prompt {prompt_by >= 0.0 && pt <= prompt_by};
		if (prompt) add(&d.prompt_charge[charge], weight);
		return prompt;
	}

	/**
	* Adds an atom absorbed after moving into a wall cell (chi = 1); see
	* Deposition::add_into_wall_cell. It goes on the element bordering
	* that wall cell whose plasma cell it came from, or the first, or is
	* unassigned if no element borders it.
	*
	* @param d The tallies.
	* @param zone Zone of the wall cell.
	* @param i x index of the wall cell.
	* @param j y index of the wall cell.
	* @param k z index of the wall cell.
	* @param prev_zone Zone of the cell the atom was in before.
	* @param prev_i x index of that cell.
	* @param prev_j y index of that cell.
	* @param pt Time since the atom was born (s).
	* @param q Its charge.
	* @param weight Its weight (atoms/s).
	* @param prompt_by Prompt window (Slots::prompt_by), below 0 for none.
	* @param v Its velocity, for the impact tallies, or nullptr.
	* @param mass Its mass (kg).
	* @param add Adder.
	* @return Index of the element it was put on, or -1.
	*/
	template <typename Add>
	FLAN_HD inline long add_into_wall_cell(const DepositionView& d,
		const int zone, const int i, const int j, const int k,
		const int prev_zone, const int prev_i, const int prev_j,
		const double pt, const int q, const double weight,
		const double prompt_by, const double* v, const double mass,
		const Add& add)
	{
		const bool prompt {add_lifetime_and_charge(d, pt, q, weight,
			prompt_by, add)};

		const long c {find_key(d.wall_cell_keys, d.n_wall_cells,
			cell_key(zone, i, j, k))};
		if (c < 0)
		{
			add_impact(d, v, mass, -1, 0.0, weight, add);
			add(&d.scalars[dep_unassigned], weight);
			return -1;
		}

		const int first {d.wall_cell_start[c]};
		const int last {d.wall_cell_start[c + 1]};
		long pick {d.wall_cell_elems[first]};
		for (int m {first}; m < last; ++m)
		{
			const int e {d.wall_cell_elems[m]};
			if (d.el.zone[e] == prev_zone && d.el.i[e] == prev_i
				&& d.el.j[e] == prev_j)
			{
				pick = e;
				break;
			}
		}
		add(&d.redeposited[pick], weight);
		if (prompt) add(&d.prompt[pick], weight);

		// The staircase element's normal, at the middle of its phi range
		add_impact(d, v, mass, pick, 0.5 * (d.el.phi_min[pick]
			+ d.el.phi_max[pick]), weight, add);
		return pick;
	}

	/**
	* Adds an atom absorbed at a wall face (see
	* Deposition::add_at_wall_face): on the element at the cell of the
	* face it crossed if there is one, otherwise only counted at the wall
	* faces. z is wrapped first if the zone is periodic in z.
	*
	* @param d The tallies.
	* @param zones The background's zones.
	* @param zone Zone whose face the atom crossed.
	* @param face The face (Zones::Face).
	* @param pos Its computational position past the face, 3 values.
	* @param pt Time since the atom was born (s).
	* @param q Its charge.
	* @param weight Its weight (atoms/s).
	* @param prompt_by Prompt window, below 0 for none.
	* @param v Its velocity, for the impact tallies, or nullptr.
	* @param mass Its mass (kg).
	* @param add Adder.
	* @return Index of the element it was put on, or -1.
	*/
	template <typename Add>
	FLAN_HD inline long add_at_wall_face(const DepositionView& d,
		const Zones::ZonesView& zones, const int zone, const int face,
		const double* pos, const double pt, const int q, const double weight,
		const double prompt_by, const double* v, const double mass,
		const Add& add)
	{
		const bool prompt {add_lifetime_and_charge(d, pt, q, weight,
			prompt_by, add)};

		if (d.n_wall_faces > 0)
		{
			// Cell along the face and toroidal cell
			const bool x_face {face == Zones::xmin || face == Zones::xmax};
			const int along {x_face ? Zones::cell_index(zones, zone, 1, pos[1])
				: Zones::cell_index(zones, zone, 0, pos[0])};
			const double z {Zones::face(zones, zone, Zones::zmin) == zone
				? Zones::wrap_z(zones, zone, pos[2]) : pos[2]};
			const int k {Zones::cell_index(zones, zone, 2, z)};
			const long c {find_key(d.wall_face_keys, d.n_wall_faces,
				wall_face_key(zone, face, along, along, k))};
			if (c >= 0)
			{
				const long e {d.wall_face_elems[c]};
				add(&d.redeposited[e], weight);
				if (prompt) add(&d.prompt[e], weight);
				add_impact(d, v, mass, e, z, weight, add);
				return e;
			}
		}

		add_impact(d, v, mass, -1, 0.0, weight, add);
		add(&d.scalars[dep_wall_face], weight);
		return -1;
	}

	/**
	* A source copied to a GPU (Source::to_device): the view's arrays are
	* on that device. Freed with free_source.
	*/
	struct SourceDevice
	{
		int device_id {-1};		///< GPU the arrays are on, -1 for none
		SourceView view {};		///< The source, with device pointers
		double weight {};		///< Atoms/s each particle stands for,
								///< source_rate / imp_num
	};

	/**
	* Deposition tallies on a GPU (Deposition::to_device): the lookups
	* copied and the tallies zeroed, on that device. Added back to the host
	* with Deposition::add_from_device and freed with free_deposition.
	*/
	struct DepositionDevice
	{
		int device_id {-1};		///< GPU the arrays are on, -1 for none
		DepositionView view {};	///< The tallies, with device pointers
	};

}  // namespace Sputtering
