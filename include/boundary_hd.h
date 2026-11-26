#pragma once
/**
* @file boundary_hd.h
*
* @brief Per-particle boundary conditions, shared by the CPU and the GPU
*        (see host_device.h).
*
* Boundary::check_bounds_cpu (boundary.cpp) and check_bounds_kernel
* (cuda/boundary.cu) apply the boundary conditions with these functions,
* in the same order: the time boundary (time_bc), the faces of the
* particle's zone (check_faces, which moves it into other zones, wraps it
* at periodic faces, puts it back at the core and absorbs it at absorbing,
* wall and lcfs faces), the wall mask chi (in_wall_cell), and the turn of
* the Cartesian velocity of a particle that went around the torus
* (turn_velocity). The CPU adds the wall deposition and self-sputtering in
* between, which the GPU doesn't have yet.
*
* Nothing here includes MPI or the standard containers, so it can be used
* in .cu files.
*/

#include <cmath>

#include "host_device.h"
#include "pcg32.h"
#include "zones_hd.h"


namespace Boundary
{
	/**
	* What happened to particles at the boundaries over a run, counted by
	* check_bounds_cpu, or on the GPU by check_bounds_kernel and read back
	* with add_fates_from_device. Reaching imp_max_time is counted
	* separately.
	*/
	struct Fates
	{
		long long crossings {};	///< Crossings into another zone
		long long absorbed {};	///< Absorbed at an absorbing face
		long long wall {};		///< Absorbed at a wall face
		long long chi {};		///< Absorbed inside the wall mask (chi)
		long long core {};		///< Put back in at a core face
		long long time_end {};	///< Reached the end of the background time
		long long runaway {};	///< Absorbed after too many crossings
		long long self_sputtered {};	///< Landed on the wall and carried
										///< on as a self-sputtered atom
		long long max_time {};	///< Followed longer than imp_max_time (set
								///< by follow_impurities, not here)

		Fates& operator+=(const Fates& other);
	};

	/// Far more zone crossings in one step than should ever happen. A
	/// particle that makes more is absorbed (see check_faces).
	inline constexpr int max_crossings {8};

	/**
	* What can happen to a particle at the boundaries in one step, as
	* indices into FateTally::n and the GPU's counters (the fields of
	* Fates, without max_time, which is counted while stepping).
	* self_sputtered is counted here only by the GPU; the CPU counts it in
	* check_bounds_cpu.
	*/
	enum FateIndex : int
	{
		fate_crossings,		///< Crossings into another zone
		fate_absorbed,		///< Absorbed at an absorbing or lcfs face
		fate_wall,			///< Absorbed at a wall face
		fate_chi,			///< Absorbed inside the wall mask (chi)
		fate_core,			///< Put back in at a core face
		fate_time_end,		///< Reached the end of the background time
		fate_runaway,		///< Absorbed after too many crossings
		fate_self_sputtered,	///< Carried on as a self-sputtered atom
		num_fate_indices
	};

	/// Counts of what happened to one particle in one step
	struct FateTally
	{
		int n[num_fate_indices] {};
	};

	/**
	* The settings of the boundary conditions that aren't in the zones:
	* the time boundary, the buffers, the lcfs, the handedness for turning
	* velocities and the wall mask.
	*/
	struct BoundaryParams
	{
		int tbound_type;	///< 0 = absorbing, 1 = periodic
		double t_min;		///< First background time (s)
		double t_max;		///< Last background time (s)
		double buffers[3];	///< imp_xbound_buffer, imp_ybound_buffer and
							///< imp_zbound_buffer
		double lcfs_x;		///< x of the last closed flux surface
		double phi_sign;	///< +1 if (R, z, Z) is right handed, -1 if not

		/**
		* The wall mask chi, flat data laid out like the zones' 3D fields
		* (Zones::flat_index with dim1 = 1), or nullptr if particles
		* aren't absorbed by it (no mask, or s3x_chi_absorb off)
		*/
		const double* chi;
	};

	/**
	* Checks if a value is outside an absorbing boundary.
	*
	* max_bound = true means bound is a maximum and the check is
	* a + buffer >= bound. max_bound = false means bound is a minimum and
	* the check is a - buffer <= bound. The buffer makes the boundary
	* trigger sooner, since the plasma can behave oddly right at the
	* boundaries. Written without branches so it is SIMD friendly.
	*
	* @param a Value to check (e.g., the particle's x coordinate).
	* @param bound Boundary value.
	* @param max_bound If true, bound is a maximum. Otherwise it is a
	*        minimum.
	* @param buffer Distance the boundary is moved inward. Defaults to 0.
	* @return 1 if the particle is outside the boundary (dead), else 0.
	*/
	FLAN_HD inline int absorbing_bc(const double a, const double bound,
		const bool max_bound, const double buffer = 0.0)
	{
		const bool ge {(a + buffer) >= bound};  // for max-bound
		const bool le {(a - buffer) <= bound};  // for min-bound
		const bool dead {max_bound ? ge : le};
		return static_cast<int>(dead);
	}

	/**
	* Applies a periodic boundary condition.
	*
	* A periodic boundary never kills a particle. A value above amax has
	* the domain length (amax - amin) subtracted, and a value below amin has
	* it added. Only a single wrap is done, so a particle is assumed to
	* move less than one domain length per time step. Written without
	* branches so it is SIMD friendly.
	*
	* @param a Value to check (e.g., the particle's y coordinate).
	* @param amin Minimum of the domain.
	* @param amax Maximum of the domain.
	* @return The wrapped value, or a unchanged if it is inside the domain.
	*/
	FLAN_HD inline double periodic_bc(const double a, const double amin,
		const double amax)
	{
		const double L {amax - amin};
		const double gt = (a > amax);   // wrap high, subtract L
		const double lt = (a < amin);   // wrap low, add L
		return a - gt * L + lt * L;
	}

	/**
	* Applies a core boundary condition.
	*
	* The coordinate being checked is a, and b and c are the other two
	* coordinates. If a is past the boundary (moved inward by a_buffer), a
	* is put back on the buffered boundary, and b and c are moved to
	* uniformly random values between their buffered minimums and
	* maximums (b drawn first). This is like the particle entering the core
	* and popping back out somewhere else along the core boundary. A
	* residence time could be added if it was useful.
	*
	* @param[in,out] a Coordinate checked against the boundary.
	* @param[in,out] b First of the other two coordinates.
	* @param[in,out] c Second of the other two coordinates.
	* @param a_min Minimum of a.
	* @param a_max Maximum of a.
	* @param a_buffer Distance the a boundary is moved inward.
	* @param b_min Minimum of b.
	* @param b_max Maximum of b.
	* @param b_buffer Buffer on each side of the b range.
	* @param c_min Minimum of c.
	* @param c_max Maximum of c.
	* @param c_buffer Buffer on each side of the c range.
	* @param max_bound If true, check the a_max side. Otherwise check the
	*        a_min side.
	* @param rng Random number generator for this particle (GPU) or
	*        thread (CPU).
	*/
	FLAN_HD inline void core_bc(double& a, double& b, double& c,
		const double a_min, const double a_max, const double a_buffer,
		const double b_min, const double b_max, const double b_buffer,
		const double c_min, const double c_max, const double c_buffer,
		const bool max_bound, pcg32& rng)
	{
		const bool hit {max_bound ? (a + a_buffer) >= a_max
			: (a - a_buffer) <= a_min};
		if (hit)
		{
			a = max_bound ? a_max - a_buffer : a_min + a_buffer;
			const double b_lo {b_min + b_buffer};
			const double b_hi {b_max - b_buffer};
			b = b_lo + rng.next_double() * (b_hi - b_lo);
			const double c_lo {c_min + c_buffer};
			const double c_hi {c_max - c_buffer};
			c = c_lo + rng.next_double() * (c_hi - c_lo);
		}
	}

	/**
	* Returns the twist-shift y offset, Delta y(x), for a particle on a
	* closed field line that wraps across the z boundary.
	*
	* This is a hook for a future twist-shift implementation and currently
	* always returns 0, so lcfs_bc does a plain periodic wrap in z. With a
	* twist-shift, a closed field line leaving z_max re-enters at z_min at
	* a different y, so Delta y(x) should match the shift Gkeyll used for
	* the background (e.g., tabulated on the x grid from a user-supplied
	* function and linearly interpolated here). The sign convention must
	* also match Gkeyll's; see lcfs_bc for how the shift is applied.
	*
	* @param x Particle's x coordinate.
	* @return Shift in y applied when wrapping from z_max to z_min (the
	*         opposite sign is applied from z_min to z_max).
	*/
	FLAN_HD inline double twist_shift_dy([[maybe_unused]] const double x)
	{
		// TODO: Implement twist-shift. Returning 0 gives a plain periodic
		// boundary in z.
		return 0.0;
	}

	/**
	* Applies a last closed flux surface (LCFS) boundary condition at one
	* z boundary.
	*
	* Inside the LCFS (x < lcfs_x) field lines are closed, so the boundary
	* is periodic: a particle past this z boundary is wrapped to the other
	* side, and its y is shifted by y_shift (the twist-shift) and wrapped
	* periodically in y. Outside the LCFS (x >= lcfs_x, the SOL) field
	* lines end on a wall, so the boundary is absorbing, including the z
	* buffer. A particle exactly at lcfs_x counts as SOL.
	*
	* ASSUMPTION: x increases outward (x < lcfs_x is the core side). This
	* is true for the geometries currently used, but may not always be.
	* If x ever increases inward, the core/SOL test here must be flipped.
	*
	* Only the side given by max_bound is checked. Crossing z_max adds
	* y_shift to y, and crossing z_min subtracts it, so a particle that
	* goes out and back returns to its original field line. The periodic
	* wrap in z uses no buffer, the same as periodic_bc.
	*
	* @param[in,out] y Particle's y coordinate. Only changed by the
	*        twist-shift on the core side.
	* @param[in,out] z Particle's z coordinate. Wrapped on the core side.
	* @param x Particle's x coordinate.
	* @param lcfs_x x coordinate of the LCFS.
	* @param y_min Minimum of y.
	* @param y_max Maximum of y.
	* @param z_min Minimum of z.
	* @param z_max Maximum of z.
	* @param z_buffer Distance the absorbing (SOL) z boundary is moved
	*        inward.
	* @param y_shift Twist-shift y offset, from twist_shift_dy. 0 gives a
	*        plain periodic boundary.
	* @param max_bound If true, check the z_max side. Otherwise check the
	*        z_min side.
	* @return 1 if the particle was absorbed (dead), else 0.
	*/
	FLAN_HD inline int lcfs_bc(double& y, double& z, const double x,
		const double lcfs_x, const double y_min, const double y_max,
		const double z_min, const double z_max, const double z_buffer,
		const double y_shift, const bool max_bound)
	{
		// Core side (closed field lines). Assumes x increases outward, see
		// above.
		if (x < lcfs_x)
		{
			const bool crossed {max_bound ? z > z_max : z < z_min};
			if (crossed)
			{
				const double L {z_max - z_min};
				z = max_bound ? z - L : z + L;

				// Twist-shift. y is wrapped since the shift can push it
				// out of the domain (y was already checked this step). With
				// y_shift = 0 this leaves an in-domain y unchanged.
				const double y_new {max_bound ? y + y_shift : y - y_shift};
				y = periodic_bc(y_new, y_min, y_max);
			}
			return 0;
		}

		// SOL side (open field lines): absorbing
		return absorbing_bc(z, max_bound ? z_max : z_min, max_bound,
			z_buffer);
	}

	/// Lower edge of a zone in a dimension
	FLAN_HD inline double zone_min(const Zones::ZonesView& zones,
		const int zone, const int dim)
	{
		return Zones::zone_edges(zones, zone, dim)[0];
	}

	/// Upper edge of a zone in a dimension
	FLAN_HD inline double zone_max(const Zones::ZonesView& zones,
		const int zone, const int dim)
	{
		return Zones::zone_edges(zones, zone, dim)[zones.n[dim][zone]];
	}

	/**
	* Puts a particle that hit a core face back in on a core face.
	*
	* The particle is put back on a core face of the same kind as the one
	* it hit (e.g., any xmin core face). If there is only one such face,
	* which is always the case for a single-zone background, core_bc puts
	* it back on that face. If there are several, as in a SOLEDGE3X core
	* split across zones, one is picked at random, weighted by its number
	* of cells, and the particle is moved to that zone before core_bc puts
	* it on the face.
	*
	* @param zones The zones.
	* @param[in,out] pos Particle position (x, y, z).
	* @param[in,out] zone Particle's zone.
	* @param face Core face the particle is being checked against.
	* @param buffers Boundary buffers in x, y and z.
	* @param rng Random number generator.
	* @return True if the particle hit the core face and was put back.
	*/
	FLAN_HD inline bool core_face(const Zones::ZonesView& zones,
		double* pos, int& zone, const int face, const double* buffers,
		pcg32& rng)
	{
		// a is the dimension normal to the face, b and c the other two
		const int a {face / 2};
		const int b {(a + 1) % 3};
		const int c {(a + 2) % 3};
		const bool max_bound {(face & 1) == 1};

		// Same check core_bc makes
		const bool hit {max_bound
			? (pos[a] + buffers[a]) >= zone_max(zones, zone, a)
			: (pos[a] - buffers[a]) <= zone_min(zones, zone, a)};

		// Pick one of several core faces, but only if the particle hit
		// this one. core_bc does the same check again below, but by then
		// the particle is on the chosen face so it always passes.
		const int* core {zones.core_zones + zones.core_offset[face]};
		const int count {zones.core_offset[face + 1]
			- zones.core_offset[face]};
		if (count > 1)
		{
			if (!hit) return false;

			// Weighted by the number of cells on the face: the first
			// whose running total is above r
			double total {};
			for (int k {}; k < count; ++k)
			{
				total += static_cast<double>(zones.n[b][core[k]])
					* zones.n[c][core[k]];
			}
			const double r {rng.next_double() * total};
			int pick {};
			double sum {};
			for (; pick < count - 1; ++pick)
			{
				sum += static_cast<double>(zones.n[b][core[pick]])
					* zones.n[c][core[pick]];
				if (r < sum) break;
			}

			zone = core[pick];
			pos[a] = max_bound ? zone_max(zones, zone, a)
				: zone_min(zones, zone, a);
		}

		core_bc(pos[a], pos[b], pos[c],
			zone_min(zones, zone, a), zone_max(zones, zone, a), buffers[a],
			zone_min(zones, zone, b), zone_max(zones, zone, b), buffers[b],
			zone_min(zones, zone, c), zone_max(zones, zone, c), buffers[c],
			max_bound, rng);
		return hit;
	}

	/**
	* Applies the time boundary: absorbing (tbound_type 0, the particle is
	* dead past the last background time) or periodic (1, wrapped between
	* the first and last), with no buffer.
	*
	* @param[in,out] t Particle's time.
	* @param p Boundary settings.
	* @param[in,out] tally Counts, added to.
	* @return 1 if the particle is dead, else 0.
	*/
	FLAN_HD inline int time_bc(double& t, const BoundaryParams& p,
		FateTally& tally)
	{
		if (p.tbound_type == 0)
		{
			const int state {absorbing_bc(t, p.t_max, true)};
			tally.n[fate_time_end] += state;
			return state;
		}
		if (p.tbound_type == 1) t = periodic_bc(t, p.t_min, p.t_max);
		return 0;
	}

	/**
	* Checks a particle against the faces of its zone, in the order xmin,
	* xmax, ymin, ymax, zmin, zmax, and applies what each face's code says
	* (see Zones::Zones):
	*
	* - face_absorbing or face_wall: absorbing_bc with the buffer of that
	*   dimension. The particle is dead and the rest is skipped.
	* - face_core: the particle is put back on a core face with core_face.
	* - face_lcfs (z faces of a single-zone background only): lcfs_bc.
	* - The zone itself (periodic): periodic_bc. Wrapping in z around a
	*   toroidal wedge moves the particle a wedge's width around the
	*   torus, which is added to turn.
	* - Another zone: if the particle is past the face, it moves into that
	*   zone. Its coordinate normal to the face is shifted from this zone's
	*   frame to the other's, and the other two are unchanged since the
	*   zones share the face's cell edges. There is no buffer. Across a z
	*   face between two toroidal zones the change in z is added to turn.
	*   The faces of the new zone are then checked from the start, so a
	*   particle can cross more than one face in a step (e.g., near an
	*   X-point). More than max_crossings crossings in a step means the
	*   zones aren't connected the way they should be or the time step is
	*   far too big, and the particle is absorbed (fate_runaway).
	*
	* @param zones The zones.
	* @param p Boundary settings.
	* @param[in,out] pos Particle position (x, y, z).
	* @param[in,out] zone Particle's zone.
	* @param[in,out] turn Toroidal angle the particle was moved by, which
	*        its velocity is turned by at the end (turn_velocity).
	* @param rng Random number generator, for the core faces.
	* @param[in,out] tally Counts, added to.
	* @param[out] wall_face The face it was absorbed at if that's a wall
	*        face, else -1.
	* @return 1 if the particle is dead, else 0.
	*/
	FLAN_HD inline int check_faces(const Zones::ZonesView& zones,
		const BoundaryParams& p, double* pos, int& zone, double& turn,
		pcg32& rng, FateTally& tally, int& wall_face)
	{
		wall_face = -1;
		int crossings {};
		int f {};
		while (f < Zones::num_faces)
		{
			const int code {Zones::face(zones, zone, f)};
			const int a {f / 2};
			const bool max_bound {(f & 1) == 1};
			const double lo {zone_min(zones, zone, a)};
			const double hi {zone_max(zones, zone, a)};

			if (code == Zones::face_absorbing || code == Zones::face_wall)
			{
				if (absorbing_bc(pos[a], max_bound ? hi : lo, max_bound,
					p.buffers[a]))
				{
					if (code == Zones::face_wall)
					{
						++tally.n[fate_wall];
						wall_face = f;
					}
					else ++tally.n[fate_absorbed];
					return 1;
				}
			}

			else if (code == Zones::face_core)
			{
				tally.n[fate_core] += core_face(zones, pos, zone, f,
					p.buffers, rng);
			}

			// LCFS: periodic inside lcfs_x, absorbing outside. Only z faces
			// can be lcfs (checked when the zones are made), so a = 2 here
			// and lcfs_bc gets y and z.
			else if (code == Zones::face_lcfs)
			{
				if (lcfs_bc(pos[1], pos[2], pos[0], p.lcfs_x,
					zone_min(zones, zone, 1), zone_max(zones, zone, 1), lo,
					hi, p.buffers[a], twist_shift_dy(pos[0]), max_bound))
				{
					++tally.n[fate_absorbed];
					return 1;
				}
			}

			// Periodic. Wrapping in z around a toroidal wedge moves the
			// particle a wedge's width around the torus.
			else if (code == zone)
			{
				const double before {pos[a]};
				pos[a] = periodic_bc(pos[a], lo, hi);
				if (a == 2 && (zones.flags[zone] & Zones::flag_wedge))
					turn += pos[a] - before;
			}

			// Another zone, entered if past the face
			else if (max_bound ? pos[a] > hi : pos[a] < lo)
			{
				if (++crossings > max_crossings)
				{
					++tally.n[fate_runaway];
					return 1;
				}

				++tally.n[fate_crossings];
				const double before {pos[a]};
				pos[a] = max_bound
					? pos[a] - hi + zone_min(zones, code, a)
					: pos[a] - lo + zone_max(zones, code, a);
				if (a == 2 && (zones.flags[zone] & Zones::flag_toroidal)
					&& (zones.flags[code] & Zones::flag_toroidal))
				{
					turn += pos[a] - before;
				}
				zone = code;
				f = 0;
				continue;
			}

			++f;
		}
		return 0;
	}

	/**
	* Checks if a particle is in a cell inside the wall (chi >= 0.5, from
	* the SOLEDGE3X wall mask).
	*
	* @param zones The zones.
	* @param chi The wall mask (BoundaryParams::chi), not nullptr.
	* @param zone Particle's zone.
	* @param pos Particle position (x, y, z).
	* @param[out] idx The cell the particle is in (x, y and z indices).
	* @return True if that cell is inside the wall.
	*/
	FLAN_HD inline bool in_wall_cell(const Zones::ZonesView& zones,
		const double* chi, const int zone, const double* pos, int* idx)
	{
		for (int d {}; d < 3; ++d) idx[d] = Zones::cell_index(zones, zone, d,
			pos[d]);
		return chi[Zones::flat_index(zones, 1, zone, 0, idx[0], idx[1],
			idx[2])] >= 0.5;
	}

	/**
	* Turns a Cartesian velocity about Z with a particle that was moved
	* around the torus by turn (radians of toroidal angle).
	*
	* @param[in,out] vX X component.
	* @param[in,out] vY Y component.
	* @param phi_sign +1 if (R, z, Z) is right handed, -1 if not.
	* @param turn Change in toroidal angle.
	*/
	FLAN_HD inline void turn_velocity(double& vX, double& vY,
		const double phi_sign, const double turn)
	{
		const double ang {phi_sign * turn};
		const double cs {std::cos(ang)};
		const double sn {std::sin(ang)};
		const double x {vX};
		vX = cs * x - sn * vY;
		vY = sn * x + cs * vY;
	}

}  // namespace Boundary
