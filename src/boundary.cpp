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
* edges. The per-particle routines are in boundary_hd.h, shared with the
* GPU version in cuda/boundary.cu.
*/
#include <array>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <vector>

#include "mpi.h"

#include "background.h"
#include "boundary.h"
#include "errors.h"
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
		self_sputtered += other.self_sputtered;
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
		std::array<long long, 9> local {fates.crossings, fates.absorbed,
			fates.wall, fates.chi, fates.core, fates.time_end,
			fates.runaway, fates.max_time, fates.self_sputtered};
		std::array<long long, 9> total {};
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
		if (total[8])
		{
			std::cout << "  " << total[8] << " of those absorbed at a wall "
				"carried on as self-sputtered atoms\n";
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
		sum.self_sputtered = total[8];
		return sum;
	}

	/**
	* The settings of the boundary conditions for check_faces and the rest
	* of boundary_hd.h, the same for the CPU and the GPU.
	*
	* @param bkg Background plasma, for the time bounds, the handedness and
	*        whether there is a wall mask.
	* @param opts Simulation options: tbound_type, the imp_*bound_buffer
	*        options, lcfs_x and s3x_chi_absorb.
	* @param chi The background's wall mask data, where the caller will use
	*        it (on the host, or bkg_d.chi on a GPU). It's only kept if the
	*        background has a mask and s3x_chi_absorb is on.
	* @return The settings.
	*/
	BoundaryParams make_params(const Background::Background& bkg,
		const Options::Options& opts, const double* chi)
	{
		const bool chi_absorb {!bkg.get_chi().get_data().empty()
			&& opts.s3x_chi_absorb_int() == 1};
		return {opts.tbound_type_int(), bkg.get_t_min(), bkg.get_t_max(),
			{opts.imp_xbound_buffer(), opts.imp_ybound_buffer(),
			opts.imp_zbound_buffer()}, opts.lcfs_x(), bkg.get_phi_sign(),
			chi_absorb ? chi : nullptr};
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
		const Zones::ZonesView& view {zones.view()};

		// Everything the shared routines (boundary_hd.h) need besides the
		// zones, the same as the GPU gets
		const BoundaryParams params {make_params(bkg, opts,
			bkg.get_chi().get_data().data())};

		// Self-sputtering, if on: particles landing on a wall element can
		// carry on as atoms they sputtered
		Sputtering::SelfSource* self {deposition != nullptr
			? deposition->self() : nullptr};

		// Set if a particle reached sput_self_max_gen (an avalanche)
		bool avalanche {false};

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

			// What happens to this particle this step
			FateTally tally {};

			// Time boundary
			double t {slots.t()[i]};
			int state {time_bc(t, params, tally)};
			slots.set_t(i, t);
			if (state)
			{
				slots.set_state(i, state);
				++local.time_end;
				continue;
			}

			// Work on local copies, written back at the end
			double pos[3] {slots.x()[i], slots.y()[i], slots.z()[i]};
			int zone {slots.zone()[i]};

			// How far the particle was moved in toroidal angle by wrapping
			// across z faces, which its velocity is turned by at the end
			double turn {};

			// Self-sputtering off element e, which the particle just landed
			// on at (pos[0], pos[1], z_wall). If it carries on as an atom it
			// sputtered, the slot gets that atom (its time goes on) and
			// true is returned.
			auto self_sputter = [&](const long e,
				const Sputtering::Impact& impact, const double z_wall)
			{
				if (self == nullptr || e < 0) return false;
				const Sputtering::Rebirth r {self->rebirth(
					static_cast<std::size_t>(e), slots.tidx()[i],
					{pos[0], pos[1], z_wall}, impact, slots.q()[i],
					slots.weight()[i], slots.generation()[i], rng)};
				if (r.avalanche)
				{
					#pragma omp atomic write
					avalanche = true;
				}
				if (!r.alive) return false;

				pos[0] = r.x;
				pos[1] = r.y;
				pos[2] = r.z;
				zone = r.zone;
				turn = 0.0;		// the new velocity is at the new z
				slots.set_xidx(i, r.i);
				slots.set_yidx(i, r.j);
				slots.set_zidx(i, r.k);
				slots.set_vX(i, r.vX);
				slots.set_vY(i, r.vY);
				slots.set_vZ(i, r.vZ);
				slots.set_vx(i, 0.0);	// set by the next Boris update
				slots.set_vy(i, 0.0);
				slots.set_vz(i, 0.0);
				slots.set_q(i, 0);
				slots.set_weight(i, r.weight);
				slots.set_pt(i, 0.0);	// time since its own sputtering
				slots.set_prompt_by(i, -1.0);
				slots.set_generation(i, slots.generation()[i] + 1);
				++local.self_sputtered;
				return true;
			};

			// The faces of its zone, and of any zones it crosses into
			int wall_face {};
			state = check_faces(view, params, pos, zone, turn, rng, tally,
				wall_face);
			if (tally.n[fate_runaway])
			{
				#pragma omp critical
				{
					std::cerr << "Warning! A particle crossed more than "
						<< max_crossings << " zone faces in one step and "
						"was absorbed. Check the zone connections or use a "
						"smaller time step.\n";
				}
			}

			// Landed on a wall face: where on it, and maybe a
			// self-sputtered atom carries on from there
			if (state && wall_face >= 0 && deposition != nullptr)
			{
				const auto face {static_cast<Zones::Face>(wall_face)};
				const Sputtering::Impact impact {{slots.vX()[i],
					slots.vY()[i], slots.vZ()[i]}, slots.mass()};
				const long e {deposition->add_at_wall_face(zones, zone, face,
					{pos[0], pos[1], pos[2]}, slots.pt()[i], slots.q()[i],
					slots.weight()[i], slots.prompt_by()[i], &impact)};
				const double z_wall {zones.get_face(zone, Zones::zmin)
					== zone ? zones.wrap_z(zone, pos[2]) : pos[2]};
				if (self_sputter(e, impact, z_wall)) state = 0;
			}

			// Absorb particles that ended up in a cell inside the wall
			int idx[3] {};
			if (params.chi != nullptr && !state
				&& in_wall_cell(view, params.chi, zone, pos, idx))
			{
				state = 1;
				++tally.n[fate_chi];
				if (deposition != nullptr)
				{
					const Sputtering::Impact impact {{slots.vX()[i],
						slots.vY()[i], slots.vZ()[i]}, slots.mass()};
					const long e {deposition->add_into_wall_cell(zone,
						idx[0], idx[1], idx[2], slots.zone()[i],
						slots.xidx()[i], slots.yidx()[i], slots.pt()[i],
						slots.q()[i], slots.weight()[i],
						slots.prompt_by()[i], &impact)};
					if (self_sputter(e, impact, pos[2])) state = 0;
				}
			}

			// The velocity is Cartesian, so it turns with the particle
			// around the torus
			if (turn != 0.0)
			{
				double vX {slots.vX()[i]};
				double vY {slots.vY()[i]};
				turn_velocity(vX, vY, params.phi_sign, turn);
				slots.set_vX(i, vX);
				slots.set_vY(i, vY);
			}

			// Write back
			slots.set_state(i, state);
			slots.set_x(i, pos[0]);
			slots.set_y(i, pos[1]);
			slots.set_z(i, pos[2]);
			slots.set_zone(i, zone);

			local.crossings += tally.n[fate_crossings];
			local.absorbed += tally.n[fate_absorbed];
			local.wall += tally.n[fate_wall];
			local.chi += tally.n[fate_chi];
			local.core += tally.n[fate_core];
			local.runaway += tally.n[fate_runaway];

		} // slot loop

		#pragma omp critical
		fates += local;

		} // omp parallel

		// A chain of self-sputtering that went on for sput_self_max_gen
		// generations: the multiplication is at or near 1, so there's no
		// steady state to find
		if (avalanche)
		{
			avalanche_error(opts.sput_self_max_gen(),
				self->chain_multiplication(), self->mean_yield(),
				self->mean_energy());
		}
	}  // check_bounds_cpu

	/**
	* Stops the run after a chain of self-sputtering went on for
	* sput_self_max_gen generations: the multiplication is at or near 1,
	* so the self-sputtered rate has no steady value to find. Called by
	* check_bounds_cpu and, for the GPU, check_bounds_wrapper.
	*
	* @param max_gen sput_self_max_gen.
	* @param chain_multiplication Multiplication along the chains so far
	*        (SelfSource::chain_multiplication).
	* @param mean_yield Mean yield of the landings so far.
	* @param mean_energy Their mean impact energy (eV).
	*/
	[[noreturn]] void avalanche_error(const int max_gen,
		const double chain_multiplication, const double mean_yield,
		const double mean_energy)
	{
		std::ostringstream msg;
		msg << "A particle was self-sputtered " << max_gen << " times "
			"(sput_self_max_gen). This looks like a self-sputtering "
			"avalanche: each atom landing on the wall sputters, on "
			"average, about as many atoms as come back, so the "
			"self-sputtered rate has no steady value. Multiplication along "
			"the chains so far: " << chain_multiplication << ", mean yield "
			<< mean_yield << " at a mean impact energy of " << mean_energy
			<< " eV.";
		Errors::fatal_error(msg.str());
	}

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
