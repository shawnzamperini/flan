/**
* @file  boundary.cpp
* @brief Routines for handling time and spatial boundary conditions.
*
* Each time step, every particle is checked against the time boundary and
* the minimum and maximum x, y and z boundaries. Each boundary can be
* absorbing (the particle is killed), periodic (the particle is wrapped to
* the other side) or, for the spatial boundaries, core (the particle is put
* back at the boundary at a random location along it). The simulation
* bounds are the first and last background times and the outer cell edges.
* The GPU version is in cuda/boundary.cu.
*/
#include <omp.h>
#include <vector>

#include "background.h"
#include "boundary.h"
#include "options.h"
#include "pcg32.h"
#include "slots.h"


namespace Boundary
{

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
	* Checks every alive particle against the time and spatial boundaries
	* and applies the boundary conditions (CPU implementation).
	*
	* The boundaries are checked in the order t, then min/max x, y and z.
	* Each boundary type is an integer from the *bound_type options:
	* 0 = absorbing, 1 = periodic and 2 = core (x, y and z only). An
	* absorbing boundary sets the slot's state and skips the remaining
	* checks as soon as the particle is dead. For a core boundary in x, y
	* and z are randomized, in y, z and x are randomized, and in z, x and y
	* are randomized. The spatial buffers are the imp_*bound_buffer
	* options, and the time boundary has no buffer. The loop over particles
	* is parallelized with OpenMP.
	*
	* @param slots Particle slots. The state, t and positions are updated
	*        in place.
	* @param bkg Background plasma, used for the simulation bounds.
	* @param opts Simulation options, used for the boundary types and
	*        buffers.
	* @param rngs Random number generators, one per OpenMP thread, used by
	*        the core boundary.
	*/
	void check_bounds_cpu(Slots::Slots& slots, 
		const Background::Background& bkg, const Options::Options& opts,
		std::vector<pcg32>& rngs)
	{

		// The different bounds are
		//   0: Absorbing
		//   1: Periodic
		//   2: Core (only x,y,z)

		#pragma omp parallel
		{

		// Grab our RNG for this thread, each thread has its own (and they're
		// seeded uniquely).
		const int tid {omp_get_thread_num()};
		pcg32& rng = rngs[tid];

		// Loop through each slot
		#pragma omp for
		for (int i = 0; i < slots.N(); i++)
		{
			// Skip dead particles
			if (slots.state()[i] > 0) continue;

			// Used in all, just define up front
			int state {};
			double new_val {};

			// --------------------
			// Time boundary
			// --------------------

			// Maximum t: Absorbing boundary
			if (opts.tbound_type_int() == 0)
			{
				state = absorbing_bc(slots.t()[i], bkg.get_t_max(), true);
				slots.set_state(i, state);

				// If dead (>0) we can skip the rest of the checks.
				if (state) continue;
			}

			// Maximum t: Periodic boundary
			else if (opts.tbound_type_int() == 1)
			{
				// We check and overwrite each time here, even though hitting
				// a BC is relatively rare. Despite this, this is still a very
				// fast process since slots.t()[i] is loaded into L1, and 
				// periodic_bc doesn't touch memory loads so the L1 cache never
				// get evicted, so it's a very fast write.
				new_val = periodic_bc(slots.t()[i], bkg.get_t_min(), 
					bkg.get_t_max());
				slots.set_t(i, new_val);
			}

			// --------------------
			// Minimum x boundary
			// --------------------
		
			// Absorbing boundary
			if (opts.min_xbound_type_int() == 0)
			{
				state = absorbing_bc(slots.x()[i], bkg.get_x_min(), false, 
					opts.imp_xbound_buffer());
				slots.set_state(i, state);

				// If dead (>0) we can skip the rest of the checks.
				if (state) continue;
			}

			// Periodic boundary
			else if (opts.min_xbound_type_int() == 1)
			{
				new_val = periodic_bc(slots.x()[i], bkg.get_x_min(), 
					bkg.get_x_max());
				slots.set_x(i, new_val);
			}

			// Core boundary
			else if (opts.min_xbound_type_int() == 2)
			{

				// Temporaries that may get updated in core_bc, which we will
				// write into slots after
				double xtmp {slots.x()[i]};
				double ytmp {slots.y()[i]};
				double ztmp {slots.z()[i]};
				core_bc(xtmp, ytmp, ztmp, bkg.get_x_min(), bkg.get_x_max(), 
					opts.imp_xbound_buffer(), bkg.get_y_min(), bkg.get_y_max(), 
					opts.imp_ybound_buffer(), bkg.get_z_min(), bkg.get_z_max(), 
					opts.imp_zbound_buffer(), false, rng);

				// Write back
				slots.set_x(i, xtmp);
				slots.set_y(i, ytmp);
				slots.set_z(i, ztmp);
			}

			// --------------------
			// Maximum x boundary
			// --------------------

			// Absorbing boundary
			if (opts.max_xbound_type_int() == 0)
			{
				state = absorbing_bc(slots.x()[i], bkg.get_x_max(), true, 
					opts.imp_xbound_buffer());
				slots.set_state(i, state);

				// If dead (>0) we can skip the rest of the checks.
				if (state) continue;
			}

			// Periodic boundary
			else if (opts.max_xbound_type_int() == 1)
			{
				new_val = periodic_bc(slots.x()[i], bkg.get_x_min(), 
					bkg.get_x_max());
				slots.set_x(i, new_val);
			}

			// Core boundary
			else if (opts.max_xbound_type_int() == 2)
			{

				// Temporaries that may get updated in core_bc, which we will
				// write into slots after
				double xtmp {slots.x()[i]};
				double ytmp {slots.y()[i]};
				double ztmp {slots.z()[i]};
				core_bc(xtmp, ytmp, ztmp, bkg.get_x_min(), bkg.get_x_max(), 
					opts.imp_xbound_buffer(), bkg.get_y_min(), bkg.get_y_max(), 
					opts.imp_ybound_buffer(), bkg.get_z_min(), bkg.get_z_max(), 
					opts.imp_zbound_buffer(), true, rng);

				// Write back
				slots.set_x(i, xtmp);
				slots.set_y(i, ytmp);
				slots.set_z(i, ztmp);
			}

			// --------------------
			// Minimum y boundary
			// --------------------
		
			// Absorbing boundary
			if (opts.min_ybound_type_int() == 0)
			{
				state = absorbing_bc(slots.y()[i], bkg.get_y_min(), false, 
					opts.imp_ybound_buffer());
				slots.set_state(i, state);

				// If dead (>0) we can skip the rest of the checks.
				if (state) continue;
			}

			// Periodic boundary
			else if (opts.min_ybound_type_int() == 1)
			{
				new_val = periodic_bc(slots.y()[i], bkg.get_y_min(), 
					bkg.get_y_max());
				slots.set_y(i, new_val);
			}

			// Core boundary
			else if (opts.min_ybound_type_int() == 2)
			{

				// Temporaries that may get updated in core_bc, which we will
				// write into slots after
				double xtmp {slots.x()[i]};
				double ytmp {slots.y()[i]};
				double ztmp {slots.z()[i]};
				core_bc(ytmp, ztmp, xtmp, bkg.get_y_min(), bkg.get_y_max(), 
					opts.imp_ybound_buffer(), bkg.get_z_min(), bkg.get_z_max(), 
					opts.imp_zbound_buffer(), bkg.get_x_min(), bkg.get_x_max(), 
					opts.imp_xbound_buffer(), false, rng);

				// Write back
				slots.set_x(i, xtmp);
				slots.set_y(i, ytmp);
				slots.set_z(i, ztmp);
			}

			// --------------------
			// Maximum y boundary
			// --------------------

			// Absorbing boundary
			if (opts.max_ybound_type_int() == 0)
			{
				state = absorbing_bc(slots.y()[i], bkg.get_y_max(), true, 
					opts.imp_ybound_buffer());
				slots.set_state(i, state);

				// If dead (>0) we can skip the rest of the checks.
				if (state) continue;
			}

			// Periodic boundary
			else if (opts.max_ybound_type_int() == 1)
			{
				new_val = periodic_bc(slots.y()[i], bkg.get_y_min(), 
					bkg.get_y_max());
				slots.set_y(i, new_val);
			}

			// Core boundary
			else if (opts.max_ybound_type_int() == 2)
			{

				// Temporaries that may get updated in core_bc, which we will
				// write into slots after
				double xtmp {slots.x()[i]};
				double ytmp {slots.y()[i]};
				double ztmp {slots.z()[i]};
				core_bc(ytmp, ztmp, xtmp, bkg.get_y_min(), bkg.get_y_max(), 
					opts.imp_ybound_buffer(), bkg.get_z_min(), bkg.get_z_max(), 
					opts.imp_zbound_buffer(), bkg.get_x_min(), bkg.get_x_max(), 
					opts.imp_xbound_buffer(), true, rng);

				// Write back
				slots.set_x(i, xtmp);
				slots.set_y(i, ytmp);
				slots.set_z(i, ztmp);
			}


			// --------------------
			// Minimum z boundary
			// --------------------
		
			// Absorbing boundary
			if (opts.min_zbound_type_int() == 0)
			{
				state = absorbing_bc(slots.z()[i], bkg.get_z_min(), false, 
					opts.imp_zbound_buffer());
				slots.set_state(i, state);

				// If dead (>0) we can skip the rest of the checks.
				if (state) continue;
			}

			// Periodic boundary
			else if (opts.min_zbound_type_int() == 1)
			{
				new_val = periodic_bc(slots.z()[i], bkg.get_z_min(), 
					bkg.get_z_max());
				slots.set_z(i, new_val);
			}

			// Core boundary
			else if (opts.min_zbound_type_int() == 2)
			{

				// Temporaries that may get updated in core_bc, which we will
				// write into slots after
				double xtmp {slots.x()[i]};
				double ytmp {slots.y()[i]};
				double ztmp {slots.z()[i]};
				core_bc(ztmp, xtmp, ytmp, bkg.get_z_min(), bkg.get_z_max(), 
					opts.imp_zbound_buffer(), bkg.get_x_min(), bkg.get_x_max(), 
					opts.imp_xbound_buffer(), bkg.get_y_min(), bkg.get_y_max(), 
					opts.imp_ybound_buffer(), false, rng);

				// Write back
				slots.set_x(i, xtmp);
				slots.set_y(i, ytmp);
				slots.set_z(i, ztmp);
			}

			// --------------------
			// Maximum z boundary
			// --------------------

			// Absorbing boundary
			if (opts.max_zbound_type_int() == 0)
			{
				state = absorbing_bc(slots.z()[i], bkg.get_z_max(), true, 
					opts.imp_zbound_buffer());
				slots.set_state(i, state);

				// If dead (>0) we can skip the rest of the checks.
				if (state) continue;
			}

			// Periodic boundary
			else if (opts.max_zbound_type_int() == 1)
			{
				new_val = periodic_bc(slots.z()[i], bkg.get_z_min(), 
					bkg.get_z_max());
				slots.set_z(i, new_val);
			}

			// Core boundary
			else if (opts.max_zbound_type_int() == 2)
			{

				// Temporaries that may get updated in core_bc, which we will
				// write into slots after
				double xtmp {slots.x()[i]};
				double ytmp {slots.y()[i]};
				double ztmp {slots.z()[i]};
				core_bc(ztmp, xtmp, ytmp, bkg.get_z_min(), bkg.get_z_max(), 
					opts.imp_zbound_buffer(), bkg.get_x_min(), bkg.get_x_max(), 
					opts.imp_xbound_buffer(), bkg.get_y_min(), bkg.get_y_max(), 
					opts.imp_ybound_buffer(), true, rng);

				// Write back
				slots.set_x(i, xtmp);
				slots.set_y(i, ytmp);
				slots.set_z(i, ztmp);
			}

		} // slot loop
		} // omp parallel
	}  // check_bounds_cpu

}  // namespace Boundary
