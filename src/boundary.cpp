/**
* @file  boundary.cpp
* @brief Routines for handling time and spatial boundary conditions.
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

	// Absorbing boundary condition. Check is value is outside the bound.
	// max_bound = true indicates that bound is a maximum, and that we are
	// checking if value >= bound. Likewise, max_bound = false checks if
	// value =< bound. buffer extends the boundary so that it triggers sooner.
	// This is because sometimes the plasma can behave oddly right at the
	// boundaries, so this can circumvent it. This function is written this way 
	// to be SIMD/vectorizable friendly for the compiler (no if's).
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

	// Periodic boundary condition. This type of BC will never kill a particle,
	// just loop it around to the other side. This function returns the new
	// looped value, or just the same value if a BC was not encountered.
	double periodic_bc(double a, double amin, double amax)
	{
		const double L = amax - amin;

		const double gt = (a > amax);   // wrap high, subtract L
		const double lt = (a < amin);   // wrap low, add L

		// Return new (or same) position
		return a - gt * L + lt * L;
	}

	// Core boundary condition. The coordinate being checked against is a. b
	// and c are the other two coordinates. If a core condition is met, b is
	// moved to a random value between b_min/b_max, likewise for c. This is
	// essentially like entering the core and then popping out at a random
	// location somewhere else along the core boundary. A residence time
	// could easily be added to this if it was useful.
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

			// Particle starts at random a/b locations (accounting for buffers)
			const double b_lo = b_min + b_buffer;
			const double b_hi = b_max - b_buffer;
			b = b_lo + rng.next_double() * (b_hi - b_lo);

			const double c_lo = c_min + c_buffer;
			const double c_hi = c_max - c_buffer;
			c = c_lo + rng.next_double() * (c_hi - c_lo);
		}
	}


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

		// Grab our RNG for this thread, each thread has it own (and they're
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

}
