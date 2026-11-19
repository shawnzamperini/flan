/**
* @file  boundary.cpp
* @brief Routines for handling time and spatial boundary conditions.
*
* Each time step, every particle is checked against the time boundary and
* the faces of its zone. The time boundary can be absorbing (the particle
* is killed) or periodic (the particle is wrapped to the other side). A
* zone face can be absorbing, wall (a material surface, which absorbs like
* an absorbing face but is counted separately), core (the particle is put
* back on a core face at a random location along it), periodic, or shared
* with another zone (the particle moves into that zone). The time bounds are the first and
* last background times, and the zone faces are the zone's outer cell
* edges. The GPU version is in cuda/boundary.cu, and only handles a single
* zone.
*/
#include <array>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <vector>

#include "mpi.h"

#include "background.h"
#include "boundary.h"
#include "indices.h"
#include "options.h"
#include "pcg32.h"
#include "slots.h"
#include "sputtering.h"
#include "zones.h"


namespace Boundary
{
	/**
	* Adds another set of counts to this one.
	*
	* @param other Counts to add.
	* @return This object.
	*/
	Fates& Fates::operator+=(const Fates& other)
	{
		crossings += other.crossings;
		absorbed += other.absorbed;
		wall += other.wall;
		chi += other.chi;
		core += other.core;
		time_end += other.time_end;
		runaway += other.runaway;
		max_time += other.max_time;
		return *this;
	}

	/**
	* Sums the particle fates over every MPI rank and prints them on rank 0.
	*
	* Must be called by every rank.
	*
	* @param fates This rank's counts.
	* @return The sums over every rank on rank 0. Other ranks get their own
	*         counts back.
	*/
	Fates print_fates(const Fates& fates)
	{
		std::array<long long, 8> local {fates.crossings, fates.absorbed,
			fates.wall, fates.chi, fates.core, fates.time_end,
			fates.runaway, fates.max_time};
		std::array<long long, 8> total {};
		MPI_Reduce(local.data(), total.data(),
			static_cast<int>(local.size()), MPI_LONG_LONG, MPI_SUM, 0,
			MPI_COMM_WORLD);

		int rank {};
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank != 0) return fates;

		std::cout << "Particle fates at the boundaries:\n"
			<< "  " << total[0] << " crossings into another zone\n"
			<< "  " << total[1] << " absorbed at an absorbing face\n"
			<< "  " << total[2] << " absorbed at a wall face\n"
			<< "  " << total[3] << " absorbed inside the wall mask (chi)\n"
			<< "  " << total[4] << " put back in at the core\n"
			<< "  " << total[5] << " reached the end of the background "
				"time\n";
		if (total[6])
		{
			std::cout << "  " << total[6] << " absorbed after crossing too "
				"many zone faces in one step\n";
		}

		Fates sum {};
		sum.crossings = total[0];
		sum.absorbed = total[1];
		sum.wall = total[2];
		sum.chi = total[3];
		sum.core = total[4];
		sum.time_end = total[5];
		sum.runaway = total[6];
		sum.max_time = total[7];
		return sum;
	}

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
	int absorbing_bc(double a, double bound, bool max_bound,
		double buffer)
	{
		// Compute both comparisons
		const bool ge {(a + buffer) >= bound};  // for max-bound
		const bool le {(a - buffer) <= bound};  // for min-bound

		// Select the correct one using max_bound as a mask
		const bool dead {max_bound ? ge : le};

		// Returns 0 if alive, 1 if dead
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
	double periodic_bc(double a, double amin, double amax)
	{
		const double L = amax - amin;

		const double gt = (a > amax);   // wrap high, subtract L
		const double lt = (a < amin);   // wrap low, add L

		// Return new (or same) position
		return a - gt * L + lt * L;
	}

	/**
	* Applies a core boundary condition.
	*
	* The coordinate being checked is a, and b and c are the other two
	* coordinates. If a is past the boundary (moved inward by a_buffer), a
	* is put back on the buffered boundary, and b and c are moved to
	* uniformly random values between their buffered minimums and
	* maximums. This is like the particle entering the core and popping
	* back out somewhere else along the core boundary. A residence time
	* could be added if it was useful.
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
	* @param rng Random number generator for this thread.
	*/
	void core_bc(double& a, double& b, double& c, double a_min,
		double a_max, double a_buffer, double b_min, double b_max,
		double b_buffer, double c_min, double c_max, double c_buffer, 
		bool max_bound, pcg32& rng)
	{
		// Only check the boundary side we were called for, accounting for 
		// a (symmetric) buffer to avoid artificial edge effects.
		const bool hit = max_bound ? (a + a_buffer) >= a_max 
			: (a - a_buffer) <= a_min;

		if (hit)
		{
			a = max_bound ? a_max - a_buffer : a_min + a_buffer;

			// Particle starts at random b/c locations (accounting for buffers)
			const double b_lo = b_min + b_buffer;
			const double b_hi = b_max - b_buffer;
			b = b_lo + rng.next_double() * (b_hi - b_lo);

			const double c_lo = c_min + c_buffer;
			const double c_hi = c_max - c_buffer;
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
	double twist_shift_dy([[maybe_unused]] double x)
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
	int lcfs_bc(double& y, double& z, double x, double lcfs_x,
		double y_min, double y_max, double z_min, double z_max,
		double z_buffer, double y_shift, bool max_bound)
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

	/**
	* Puts a particle that hit a core face back in on a core face.
	*
	* The particle is put back on a core face of the same kind as the one
	* it hit (e.g., any xmin core face). If there is only one such face,
	* which is always the case for a single-zone background, core_bc puts
	* it back on that face exactly as it always has. If there are several,
	* as in a SOLEDGE3X core split across zones, one is picked at random,
	* weighted by its number of cells, and the particle is moved to that
	* zone before core_bc puts it on the face.
	*
	* @param pos Particle position (x, y, z), updated in place.
	* @param zone Particle's zone, updated in place.
	* @param face Core face the particle is being checked against.
	* @param core_zones Zones with a core face of the same kind, from
	*        Zones::Zones::get_core_zones(face).
	* @param zones The zones.
	* @param buffers Boundary buffers in x, y and z.
	* @param rng Random number generator for this thread.
	* @return True if the particle hit the core face and was put back.
	*/
	bool core_face(std::array<double, 3>& pos, int& zone,
		const Zones::Face face, const std::vector<int>& core_zones,
		const Zones::Zones& zones, const std::array<double, 3>& buffers,
		pcg32& rng)
	{
		// a is the dimension normal to the face, b and c the other two
		const int a {face / 2};
		const int b {(a + 1) % 3};
		const int c {(a + 2) % 3};
		const bool max_bound {(face & 1) == 1};

		// Same check core_bc makes
		const bool hit {max_bound
			? (pos[a] + buffers[a]) >= zones.get_max(zone, a)
			: (pos[a] - buffers[a]) <= zones.get_min(zone, a)};

		// Pick one of several core faces, but only if the particle hit
		// this one. core_bc does the same check again below, but by then
		// the particle is on the chosen face so it always passes.
		if (core_zones.size() > 1)
		{
			if (!hit) return false;

			// Weighted by the number of cells on the face
			std::vector<double> weights {};
			double total {};
			for (const int z : core_zones)
			{
				total += static_cast<double>(zones.get_n(z, b))
					* zones.get_n(z, c);
				weights.push_back(total);
			}
			const double r {rng.next_double() * total};
			std::size_t pick {};
			while (pick < weights.size() - 1 && r >= weights[pick]) ++pick;

			zone = core_zones[pick];
			pos[a] = max_bound ? zones.get_max(zone, a)
				: zones.get_min(zone, a);
		}

		core_bc(pos[a], pos[b], pos[c],
			zones.get_min(zone, a), zones.get_max(zone, a), buffers[a],
			zones.get_min(zone, b), zones.get_max(zone, b), buffers[b],
			zones.get_min(zone, c), zones.get_max(zone, c), buffers[c],
			max_bound, rng);
		return hit;
	}

	/**
	* Checks every alive particle against the time boundary and the faces
	* of its zone, and applies the boundary conditions (CPU
	* implementation).
	*
	* The time boundary comes from tbound_type: 0 = absorbing, 1 =
	* periodic, with no buffer. Then the faces of the particle's zone are
	* checked in the order xmin, xmax, ymin, ymax, zmin, zmax. What happens
	* at a face depends on its code (see Zones::Zones):
	*
	* - face_absorbing or face_wall: absorbing_bc with the
	*   imp_*bound_buffer of that dimension. A dead particle skips the
	*   remaining checks. The two are counted separately in fates.
	* - face_core: the particle is put back on a core face with core_face.
	* - face_lcfs (z faces of a single-zone background only): lcfs_bc,
	*   periodic in z (with the twist-shift in y) for particles inside
	*   lcfs_x and absorbing, with imp_zbound_buffer, outside it.
	* - The zone itself (periodic): periodic_bc. If the zone is a toroidal
	*   wedge (Zones::Zones::is_wedge), wrapping in z moves the particle a
	*   wedge's width around the torus, so its Cartesian velocity (vX, vY)
	*   is rotated about Z by phi_sign times the change in z. The same is
	*   done for a z face between two toroidal zones.
	* - Another zone: if the particle is past the face, it moves into that
	*   zone. Its coordinate normal to the face is shifted from this zone's
	*   frame to the other's, and the other two are unchanged since the
	*   zones share the face's cell edges. There is no buffer. The faces of
	*   the new zone are then checked from the start, so a particle can
	*   cross more than one face in a step (e.g., near an X-point).
	*
	* A SOLEDGE3X background has a wall mask chi. If s3x_chi_absorb is on,
	* a particle that ends up in a cell with chi >= 0.5 (inside the wall) is
	* absorbed.
	*
	* For a single-zone background whose faces come from the *bound_type
	* options, this does exactly what checking each bound_type in turn did
	* before zones. A particle that crosses more than max_crossings faces in
	* one step is absorbed with a warning, since that means the zones
	* aren't connected the way they should be or the time step is far too
	* big. The loop over particles is parallelized with OpenMP.
	*
	* @param slots Particle slots. The state, t, position and zone are
	*        updated in place.
	* @param bkg Background plasma, used for the time bounds and the zones.
	* @param opts Simulation options, used for the time boundary type and
	*        the spatial buffers.
	* @param rngs Random number generators, one per OpenMP thread, used by
	*        the core boundary.
	* @param fates Counts of what happened to particles at the boundaries,
	*        added to.
	* @param deposition If given, particles absorbed in a wall cell or at a
	*        wall face are added to it (where on the wall, how long after
	*        they were born, and their charge). At a wall face, where on it
	*        comes from the particle's position as it crosses. The cell a
	*        particle was in before this step, which it came into the wall
	*        from, is still in its slot's cell indices, since they're only
	*        updated after the boundary check.
	*/
	void check_bounds_cpu(Slots::Slots& slots, 
		const Background::Background& bkg, const Options::Options& opts,
		std::vector<pcg32>& rngs, Fates& fates,
		Sputtering::Deposition* deposition)
	{
		const Zones::Zones& zones {bkg.get_zones()};
		const std::array<double, 3> buffers {opts.imp_xbound_buffer(),
			opts.imp_ybound_buffer(), opts.imp_zbound_buffer()};

		// Zones with a core face, for each kind of face
		std::array<std::vector<int>, Zones::num_faces> core_zones {};
		for (int f {}; f < Zones::num_faces; ++f)
		{
			core_zones[f] = zones.get_core_zones(static_cast<Zones::Face>(f));
		}

		// Handedness of (R, z, Z), for turning velocities on a wedge
		const double phi_sign {bkg.get_phi_sign()};

		// Far more crossings in one step than should ever happen
		constexpr int max_crossings {8};

		// SOLEDGE3X wall mask, if particles are absorbed by the wall
		const auto& chi {bkg.get_chi()};
		const bool chi_absorb {!chi.get_data().empty()
			&& opts.s3x_chi_absorb_int() == 1};

		#pragma omp parallel
		{

		// Grab our RNG for this thread, each thread has its own (and they're
		// seeded uniquely).
		const int tid {omp_get_thread_num()};
		pcg32& rng = rngs[tid];

		// This thread's counts, added to fates at the end
		Fates local {};

		// Loop through each slot
		#pragma omp for
		for (int i = 0; i < slots.N(); i++)
		{
			// Skip dead particles
			if (slots.state()[i] > 0) continue;

			// --------------------
			// Time boundary
			// --------------------

			// Maximum t: Absorbing boundary
			if (opts.tbound_type_int() == 0)
			{
				const int state {absorbing_bc(slots.t()[i], bkg.get_t_max(), 
					true)};
				slots.set_state(i, state);

				// If dead (>0) we can skip the rest of the checks.
				if (state)
				{
					++local.time_end;
					continue;
				}
			}

			// Maximum t: Periodic boundary
			else if (opts.tbound_type_int() == 1)
			{
				// We check and overwrite each time here, even though hitting
				// a BC is relatively rare. Despite this, this is still a very
				// fast process since slots.t()[i] is loaded into L1, and 
				// periodic_bc doesn't touch memory loads so the L1 cache never
				// get evicted, so it's a very fast write.
				slots.set_t(i, periodic_bc(slots.t()[i], bkg.get_t_min(), 
					bkg.get_t_max()));
			}

			// --------------------
			// Zone faces
			// --------------------

			// Work on local copies, written back at the end
			std::array<double, 3> pos {slots.x()[i], slots.y()[i], 
				slots.z()[i]};
			int zone {slots.zone()[i]};
			int state {};
			int crossings {};

			// How far the particle was moved in toroidal angle by wrapping
			// across z faces, which its velocity is turned by at the end
			double turn {};

			// Check each face in turn. Crossing into another zone starts
			// over at the first face of the new zone.
			int f {};
			while (f < Zones::num_faces)
			{
				const auto face {static_cast<Zones::Face>(f)};
				const int code {zones.get_face(zone, face)};
				const int a {f / 2};
				const bool max_bound {(f & 1) == 1};
				const double lo {zones.get_min(zone, a)};
				const double hi {zones.get_max(zone, a)};

				if (code == Zones::face_absorbing || code == Zones::face_wall)
				{
					state = absorbing_bc(pos[a], max_bound ? hi : lo, 
						max_bound, buffers[a]);
					if (state)
					{
						++(code == Zones::face_wall ? local.wall
							: local.absorbed);
						if (deposition != nullptr && code == Zones::face_wall)
						{
							deposition->add_at_wall_face(zones, zone, face,
								pos, slots.pt()[i], slots.q()[i],
								slots.weight()[i], slots.prompt_by()[i]);
						}
						break;
					}
				}

				else if (code == Zones::face_core)
				{
					local.core += core_face(pos, zone, face, core_zones[f],
						zones, buffers, rng);
				}

				// LCFS: periodic inside lcfs_x, absorbing outside. Only z
				// faces can be lcfs (checked when the zones are made), so
				// a = 2 here and lcfs_bc gets y and z.
				else if (code == Zones::face_lcfs)
				{
					state = lcfs_bc(pos[1], pos[2], pos[0], opts.lcfs_x(),
						zones.get_min(zone, 1), zones.get_max(zone, 1), lo,
						hi, buffers[a], twist_shift_dy(pos[0]), max_bound);
					if (state)
					{
						++local.absorbed;
						break;
					}
				}

				// Periodic. Wrapping in z around a toroidal wedge moves the
				// particle a wedge's width around the torus.
				else if (code == zone)
				{
					const double before {pos[a]};
					pos[a] = periodic_bc(pos[a], lo, hi);
					if (a == 2 && zones.is_wedge(zone))
						turn += pos[a] - before;
				}

				// Another zone, entered if past the face
				else if (max_bound ? pos[a] > hi : pos[a] < lo)
				{
					if (++crossings > max_crossings)
					{
						#pragma omp critical
						{
							std::cerr << "Warning! A particle crossed more "
								"than " << max_crossings << " zone faces in "
								"one step and was absorbed. Check the zone "
								"connections or use a smaller time step.\n";
						}
						state = 1;
						++local.runaway;
						break;
					}

					++local.crossings;
					const double before {pos[a]};
					pos[a] = max_bound 
						? pos[a] - hi + zones.get_min(code, a)
						: pos[a] - lo + zones.get_max(code, a);
					if (a == 2 && zones.is_toroidal(zone)
						&& zones.is_toroidal(code))
					{
						turn += pos[a] - before;
					}
					zone = code;
					f = 0;
					continue;
				}

				++f;
			}

			// Absorb particles that ended up in a cell inside the wall
			if (chi_absorb && !state)
			{
				const int xidx {Indices::get_nearest_cell_index_cpu(
					zones.get_grid_x(zone), pos[0])};
				const int yidx {Indices::get_nearest_cell_index_cpu(
					zones.get_grid_y(zone), pos[1])};
				const int zidx {Indices::get_nearest_cell_index_cpu(
					zones.get_grid_z(zone), pos[2])};
				if (chi(zone, xidx, yidx, zidx) >= 0.5)
				{
					state = 1;
					++local.chi;
					if (deposition != nullptr)
					{
						deposition->add_into_wall_cell(zone, xidx, yidx,
							zidx, slots.zone()[i], slots.xidx()[i],
							slots.yidx()[i], slots.pt()[i], slots.q()[i],
							slots.weight()[i], slots.prompt_by()[i]);
					}
				}
			}

			// The velocity is Cartesian, so it turns with the particle
			// around the torus
			if (turn != 0.0)
			{
				const double ang {phi_sign * turn};
				const double cs {std::cos(ang)};
				const double sn {std::sin(ang)};
				const double vX {slots.vX()[i]};
				const double vY {slots.vY()[i]};
				slots.set_vX(i, cs * vX - sn * vY);
				slots.set_vY(i, sn * vX + cs * vY);
			}

			// Write back
			slots.set_state(i, state);
			slots.set_x(i, pos[0]);
			slots.set_y(i, pos[1]);
			slots.set_z(i, pos[2]);
			slots.set_zone(i, zone);

		} // slot loop

		#pragma omp critical
		fates += local;

		} // omp parallel
	}  // check_bounds_cpu

	/**
	* check_bounds_cpu without keeping the counts of what happened to the
	* particles.
	*
	* @param slots Particle slots.
	* @param bkg Background plasma.
	* @param opts Simulation options.
	* @param rngs Random number generators, one per OpenMP thread.
	*/
	void check_bounds_cpu(Slots::Slots& slots, 
		const Background::Background& bkg, const Options::Options& opts,
		std::vector<pcg32>& rngs)
	{
		Fates fates {};
		check_bounds_cpu(slots, bkg, opts, rngs, fates);
	}

}  // namespace Boundary
