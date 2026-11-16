/**
* @file  collisions.cu
*
* @brief CUDA routines for updating particle velocity due to Coulomb
*        collisions via the Nanbu algorithm.
*
* GPU version of collisions.cpp. Derivation and steps taken from:
* Nanbu, K. Theory of cumulative small-angle collisions in plasmas.
* Phys. Rev. E 55, 4642–4652 (1997).
*
* Each time step, the cumulative effect of all the small-angle collisions
* with a background species is replaced by a single deflection of the
* relative velocity by an angle chi. How collisional the step is is set by
* s, which sets the shape A(s) of the distribution chi is sampled from.
*/
#include <cuda_runtime.h>
#include <cstdio>
#include <math.h>

#include "background_device.h"
#include "constants.h"
#include "options_device.h"
#include "pcg32.h"
#include "slots_device.h"

#include "impurity_stats.cuh"
#include "interpolate.cuh"
#include "nanbu_s_a.cuh"
#include "utilities.cuh"


namespace Collisions
{
	/**
	* Values returned from nanbu_calc_s_cuda: the s term and the Cartesian
	* components of the relative velocity (m/s).
	*/
	struct NanbuSResult
	{
		double s;
		double gX;
		double gY;
		double gZ;
	};

	/**
	* Simple struct to hold the three components of a vector.
	*/
	struct Vec3d
	{
		double x;
		double y;
		double z;
	};


	/**
	* Calculates the electron Debye length,
	* lambda_D = sqrt(eps0 kTe / (ne e^2)).
	*
	* @param te Electron temperature (eV).
	* @param ne Electron density (m^-3).
	* @return Debye length (m).
	*/
	__device__ __forceinline__
	double calc_debye_length_cuda(const double te, const double ne)
	{
		return sqrt(Constants::eps0 * te * Constants::ev_to_j / 
			(ne * Constants::charge_e * Constants::charge_e)); 
	}


	/**
	* Samples the velocity of a background particle from a drifting
	* Maxwellian.
	*
	* Each component is drawn from a normal distribution with a mean equal
	* to the flow velocity and a standard deviation of sqrt(kT/m). We call
	* this thermal sampling of the background velocity.
	*
	* @param T Temperature of the background species (eV).
	* @param uX Cartesian X flow velocity of the background species (m/s).
	* @param uY Cartesian Y flow velocity of the background species (m/s).
	* @param uZ Cartesian Z flow velocity of the background species (m/s).
	* @param m Mass of the background species (kg).
	* @param rng Random number generator for this particle.
	* @param[out] bkg_vX Sampled Cartesian X velocity (m/s).
	* @param[out] bkg_vY Sampled Cartesian Y velocity (m/s).
	* @param[out] bkg_vZ Sampled Cartesian Z velocity (m/s).
	*/
	__device__ __forceinline__
	void sample_bkg_velocity_cuda(const double T, const double uX, const double uY, 
		const double uZ, const double m, pcg32& rng, double& bkg_vX, 
		double& bkg_vY, double& bkg_vZ)
	{
		// Sampling from a drifting Maxwellian, which can be done by just
		// sampling from a normally distributed number with uXYZ = mean = drift 
		// and mu = sqrt(kT/m). We call this thermal sampling of the background
		// velocity. Unfortunately throwing away the second random number 
		// here. T [eV], m [kg]
		double mu {sqrt(T * Constants::ev_to_j / m)};  // m/s
		bkg_vX = rng.normal(uX, mu);
		bkg_vY = rng.normal(uY, mu);
		bkg_vZ = rng.normal(uZ, mu);
	}


	/**
	* Device version of Collisions::nanbu_s_coeff: the part of the s term
	* in the Nanbu collision model (Eq. 19) that is the same for every
	* background particle the impurity could collide with at this location,
	* s = coeff * dt / g^3. It is calculated once per particle and time
	* step, however many sub-steps the collision is split into.
	*
	* The Coulomb logarithm uses the Debye length and the expectation
	* value of the impact parameter b0, and is held at or above 1. The
	* background species is assumed to be singly charged with ni = ne. A
	* message is printed if the Coulomb logarithm is NaN.
	*
	* @param q Charge state of the impurity. Must not be 0.
	* @param T Temperature of the background species colliding with the
	*        impurity (eV).
	* @param Te Electron temperature (eV), used for the Debye length.
	* @param ne Electron density (m^-3).
	* @param imp_mass Mass of the impurity (kg).
	* @param mass_kg Mass of the background species (kg).
	* @return coeff (m^3/s^4), so that s = coeff * dt / g^3.
	*/
	__device__ __forceinline__
	double nanbu_s_coeff_cuda(const int q, const double T, const double Te,
		const double ne, const double imp_mass, const double mass_kg)
	{
		// Reduced mass
		const double mu_ab = imp_mass * mass_kg / (imp_mass + mass_kg);

		double expect_g_sq = 3.0 * T * Constants::ev_to_j / mu_ab;

		expect_g_sq = fmax(expect_g_sq, 1.0e-3);

		// Impact parameter, b0
		const double expect_b0 = q * Constants::charge_e * Constants::charge_e /
			(2.0 * Constants::pi * Constants::eps0 * mu_ab * expect_g_sq);

		// Debye length
		const double debye_length = calc_debye_length_cuda(Te, ne);

		// Coloumb logarithm, checking for NaNs and not letting it go below 1
		double ln_alpha = log(debye_length / expect_b0);
		if (isnan(ln_alpha)) printf("nan ln_alpha\n");
		ln_alpha = fmax(ln_alpha, 1.0);

		// Everything in s but dt / g^3
		const double square_term = q * Constants::charge_e * 
			Constants::charge_e / (Constants::eps0 * mu_ab);
		return ln_alpha / (4.0 * Constants::pi) * square_term *
			square_term * ne;

	}  // nanbu_s_coeff_cuda


	/**
	* Samples the background particle an impurity collides with from a
	* drifting Maxwellian (flow plus thermal sampling) and calculates the
	* s term in the Nanbu collision model (Eq. 19) for that pair, s =
	* s_coeff * dt / g^3, with g their relative speed. The instantaneous g
	* is used for both its direction and its magnitude, as in Nanbu's
	* paper.
	*
	* Each call samples a new background particle, so calling it once per
	* sub-step is what lets the impurity meet several of them in one time
	* step.
	*
	* @param s_coeff From nanbu_s_coeff_cuda, at the impurity location.
	* @param dt Time this collision covers (s): the time step, or one
	*        sub-step of it.
	* @param T Temperature of the background species colliding with the
	*        impurity (eV).
	* @param mass_kg Mass of the background species (kg).
	* @param uX Cartesian X flow velocity of the background at the impurity
	*        location (m/s).
	* @param uY Cartesian Y flow velocity (m/s).
	* @param uZ Cartesian Z flow velocity (m/s).
	* @param vX Cartesian X velocity of the impurity (m/s).
	* @param vY Cartesian Y velocity of the impurity (m/s).
	* @param vZ Cartesian Z velocity of the impurity (m/s).
	* @param rng Random number generator for this particle.
	* @return The s term and the Cartesian components of the relative
	*         velocity (m/s).
	*/
	__device__ __forceinline__
	NanbuSResult nanbu_calc_s_cuda(const double s_coeff, const double dt,
		const double T, const double mass_kg, const double uX,
		const double uY, const double uZ, const double vX, const double vY,
		const double vZ, pcg32& rng)
	{
		// Sample background particle velocity
		double bkg_vX, bkg_vY, bkg_vZ;
		sample_bkg_velocity_cuda(T, uX, uY, uZ, mass_kg, rng, bkg_vX, bkg_vY, 
			bkg_vZ);

		// Instantaneous relative velocity
		double inst_gX = vX - bkg_vX;
		double inst_gY = vY - bkg_vY;
		double inst_gZ = vZ - bkg_vZ;

		// Magnitude
		double inst_g = sqrt(inst_gX * inst_gX + inst_gY * inst_gY +
				 inst_gZ * inst_gZ);
		inst_g = fmax(inst_g, 1.0e-3);

		// Calculate s (Eq. 19)
		const double s = s_coeff / (inst_g * inst_g * inst_g) * dt;

		// Return as struct
		return {s, inst_gX, inst_gY, inst_gZ};

	}  // nanbu_calc_s_cuda


	/**
	* Device version of Collisions::nanbu_s0: the Nanbu s term for a
	* background particle moving at its thermal speed v_th = sqrt(T / m)
	* relative to the impurity, which measures how collisional the time
	* step is.
	*
	* @param s_coeff From nanbu_s_coeff_cuda.
	* @param dt Time step (s).
	* @param T Temperature of the background species (eV).
	* @param mass_kg Mass of the background species (kg).
	* @return s0 (dimensionless).
	*/
	__device__ __forceinline__
	double nanbu_s0_cuda(const double s_coeff, const double dt,
		const double T, const double mass_kg)
	{
		const double v_th = sqrt(T * Constants::ev_to_j / mass_kg);
		return s_coeff * dt / (v_th * v_th * v_th);
	}


	/**
	* Device version of Collisions::nanbu_num_substeps (see collisions.h
	* for why collisions are sub-cycled): the number of sub-steps to split
	* one Nanbu collision into so each has s0 at or below target, at most
	* cap. 1 if s0 is at or below target or isn't a number.
	*
	* @param s0 From nanbu_s0_cuda, for the full time step.
	* @param target Largest s0 a sub-step is meant to have
	*        (imp_coll_subcycle_target).
	* @param cap Largest number of sub-steps (imp_coll_subcycle_cap).
	* @return Number of sub-steps, from 1 to cap.
	*/
	__device__ __forceinline__
	int nanbu_num_substeps_cuda(const double s0, const double target,
		const int cap)
	{
		// Also catches NaN, for which the comparison is false
		if (!(s0 > target)) return 1;

		// Compare as doubles first: s0 / target can be far too large for
		// an int in very collisional plasma
		const double n = ceil(s0 / target);
		if (n >= cap) return cap;
		const int n_int = static_cast<int>(n);
		return n_int < 1 ? 1 : n_int;
	}


	/**
	* Calculates the A term in the Nanbu collision model (Eq. 13), which sets
	* the shape of the distribution the deflection angle is sampled from.
	*
	* For s < 0.01 (weakly collisional) A = 1 / s, and for s > 5 (strongly
	* collisional) A = 3 exp(-s). In between, A is linearly interpolated
	* from the precomputed table in nanbu_s_a.cuh, which was calculated
	* using python/calc_nanbu.py.
	*
	* @param s The s term from nanbu_calc_s_cuda.
	* @return The A term.
	*/
	__device__
	double nanbu_calc_A_cuda(const double s)
	{
		// Not very collisional approximation
		if (s < 0.01)
		{
			return 1.0 / s;
		}

		// Very collisional approximation
		else if (s > 5.0)
		{
			return 3.0 * exp(-s);
		}

		// Little bit of collisional, little bit of not collisional
		else
		{
			// Use interpolation function to find A value from precomputed 
			// arrays of A(s) values found in nanbu_s_a.cuhh. These values were 
			// calculated using python/calc_nanbu.py.
			// Like the CPU code, we use a uniform spacing interpolate 
			// because it's way faster. 
			//return Utilities::linear_interpolate_cuda(Nanbu::s, Nanbu::A, s);
			return Utilities::linear_interpolate_uniform_cuda(Nanbu::s, 
				Nanbu::A, s);
		}
	}  // nanbu_calc_A


	/**
	* Samples the deflection angle chi in the Nanbu collision model (Eq. 17).
	*
	* cos(chi) is found by directly inverting the cumulative distribution
	* built from A, using a uniform random number U. For s < 0.02 the small-s
	* limit 1 + s ln(U) is used, and for s > 6 the deflection is isotropic,
	* 2U - 1. cos(chi) is clamped to [-1, 1] to guard against rounding
	* errors.
	*
	* @param s The s term from nanbu_calc_s_cuda.
	* @param A The A term from nanbu_calc_A_cuda.
	* @param rng Random number generator for this particle.
	* @return Deflection angle chi (radians).
	*/
	__device__ 
	double nanbu_calc_chi_cuda(const double s, const double A, pcg32& rng)
	{
		// Uniform random number in (0,1]
		double U = rng.next_double();
		U = fmax(U, 1.0e-12);

		// Compute all candidate values
		const double cos_small = 1.0 + s * log(U);
		const double cos_large = 2.0 * U - 1.0;

		const double expA  = exp(-A);
		const double sinhA = sinh(A);

		const double cos_mid = log(expA + 2.0 * U * sinhA) / A;

		// Select branch
		double cos_chi = cos_mid;

		if (s < 0.02) cos_chi = cos_small;
		else if (s > 6.0) cos_chi = cos_large;

		// Clamp
		cos_chi = fmin(1.0, fmax(-1.0, cos_chi));

		return acos(cos_chi);
	}  // nanbu_calc_chi_cuda

	
	/**
	* Calculates the post-collision Cartesian velocity of the impurity in the
	* Nanbu collision model (Eq. 20a).
	*
	* The relative velocity g is rotated by the deflection angle chi about a
	* random azimuthal angle, uniformly distributed in [0, 2 pi). The case
	* where g is parallel to X is handled separately to avoid dividing by
	* zero.
	*
	* @param gX Cartesian X component of the relative velocity (m/s).
	* @param gY Cartesian Y component of the relative velocity (m/s).
	* @param gZ Cartesian Z component of the relative velocity (m/s).
	* @param chi Deflection angle from nanbu_calc_chi_cuda (radians).
	* @param mass_a Mass of the impurity (kg).
	* @param mass_b Mass of the background species (kg).
	* @param vX Cartesian X velocity of the impurity before the collision
	*        (m/s).
	* @param vY Cartesian Y velocity of the impurity before the collision
	*        (m/s).
	* @param vZ Cartesian Z velocity of the impurity before the collision
	*        (m/s).
	* @param rng Random number generator for this particle.
	* @return Post-collision velocity in Cartesian space (m/s).
	*/
	__device__ 
	Vec3d nanbu_post_coll_cuda(const double gX, const double gY,
		const double gZ, const double chi, const double mass_a, 
		const double mass_b, const double vX, const double vY, const double vZ,
		pcg32& rng)
	{
		// For calculating h components
		const double g_perp = sqrt(gY*gY + gZ*gZ);
		const double g = sqrt(gX*gX + gY*gY + gZ*gZ);
		const double eps =2.0 * Constants::pi * rng.next_double();
		const double cos_eps = cos(eps);
		const double sin_eps = sin(eps);

		// Calculate the h components. If g_perp = 0 then g is
		// parallel to X and can cause NaNs/infs when we divide by g_perp, so
		// we handle accordingly.
		double hX;
		double hY;
		double hZ;

		if (g_perp < Constants::small)
		{
			hX = 0.0;
			hY = g * cos_eps;
			hZ = g * sin_eps;
		}
		else
		{
			hX = g_perp * cos_eps;
			hY = -(gY * gX * cos_eps + g  * gZ * sin_eps) / g_perp;
			hZ = -(gZ * gX * cos_eps - g  * gY * sin_eps) / g_perp;
		}

		// Calculate post collision velocities
		const double mu = mass_b / (mass_a + mass_b);
		const double cos_chi = cos(chi);
		const double sin_chi = sin(chi);

		Vec3d v;

		// From Nanbu paper, assuming self-consistent particle interaction
		v.x = vX - mu * (gX * (1.0 - cos_chi) + hX * sin_chi);
		v.y = vY - mu * (gY * (1.0 - cos_chi) + hY * sin_chi);
		v.z = vZ - mu * (gZ * (1.0 - cos_chi) + hZ * sin_chi);

		// Return struct
		return v;

	}  // nanbu_post_coll_cuda


	/**
	* Kernel that updates the Cartesian velocity of each particle due to
	* Coulomb collisions with one background species using the Nanbu model.
	*
	* One thread handles one slot, and dead particles are skipped. ne, Te
	* (and Ti for ions) are interpolated at the particle location and held
	* at or above min_ne, min_te and min_ti. Then s, A and chi are calculated
	* and the velocity is updated with nanbu_post_coll_cuda. Since this
	* happens after the Boris update, the updated velocity is v at n+1/2. A
	* message is printed if chi is NaN.
	*
	* Where the time step is too collisional for the model, the collision
	* is split into sub-steps, each with a newly sampled background
	* particle, set by imp_coll_subcycle_target and imp_coll_subcycle_cap
	* (see the notes on sub-cycling in collisions.h).
	*
	* @param slots_d Particle slots on the device. The velocities are
	*        updated in place.
	* @param bkg_d Background plasma on the device.
	* @param elec If true, collide with electrons (Te). Otherwise collide
	*        with ions (Ti).
	* @param dt Time step (s).
	* @param stats_d Impurity statistics on the device. Currently unused,
	*        since recording s is commented out.
	* @param rngs_d Random number generators on the device, one per slot.
	* @param elec_mass_amu Electron mass (amu).
	* @param ion_mass_amu Main ion mass (amu).
	* @param opts_d Simulation options on the device, used for min_ne,
	*        min_te, min_ti, imp_coll_subcycle_target and
	*        imp_coll_subcycle_cap.
	*/
	__global__
	void nanbu_coll_kernel(Slots::SlotsDevice slots_d, 
		const Background::BackgroundDevice bkg_d, bool elec, const double dt, 
		ImpurityStats::StatisticsDevice stats_d, pcg32* rngs_d, 
		const double elec_mass_amu, const double ion_mass_amu, 
		const Options::OptionsDevice* opts_d)
	{
		// Derivation and steps taken from:
		// Nanbu, K. Theory of cumulative small-angle collisions in plasmas. 
		// Phys. Rev. E 55, 4642–4652 (1997).

		// Global index
		int i = blockIdx.x * blockDim.x + threadIdx.x;

		// Don't try and access beyond the number of slots (segfault)
		if (i >= slots_d.N) return;

		// If particle isn't alive then skip
		if (slots_d.state[i] > 0) return;

		// Neutral particles do not scatter, so skip. Without this you will
		// get an crash (or hang) anyways.
		int q {slots_d.q[i]};
		if (q == 0) return;

		// Local variables
		int tidx {slots_d.tidx[i]};
		int xidx {slots_d.xidx[i]};
		int yidx {slots_d.yidx[i]};
		int zidx {slots_d.zidx[i]};
		double t {slots_d.t[i]};
		double x {slots_d.x[i]};
		double y {slots_d.y[i]};
		double z {slots_d.z[i]};
		double vX {slots_d.vX[i]};
		double vY {slots_d.vY[i]};
		double vZ {slots_d.vZ[i]};

		// Interpolate for ne, Te
		// Get interpolation stencil at this location for 4D arrays
		Interpolate::InterpolationStencil stencil_4d
			{Interpolate::build_stencil_4d(bkg_d, tidx, xidx, yidx, zidx, t, x, 
			y, z)};

		// Interpolate to get background flow at particle location
		double ne {Interpolate::interpolate_field_4d(bkg_d.ne, 
			stencil_4d, t, x, y, z)};
		double Te {Interpolate::interpolate_field_4d(bkg_d.te, 
			stencil_4d, t, x, y, z)};

		// Don't let them go below zero. This can happen near the edge of grids.
		ne = fmax(ne, opts_d->min_ne);
		Te = fmax(Te, opts_d->min_te);

		// Load some reusable variables based on which species. If elec = true,
		// then electrons and ions if not. 
		double T {};
		double mass_kg {};
		if (elec)
		{
			mass_kg = elec_mass_amu * Constants::amu_to_kg;	
			T = Te;
		}
		else
		{
			mass_kg = ion_mass_amu * Constants::amu_to_kg;	
			T = Interpolate::interpolate_field_4d(bkg_d.ti, stencil_4d, t, x, 
				y, z);
			T = fmax(opts_d->min_ti, T);
		}

		// Everything in s but dt / g^3 (see nanbu_s_coeff_cuda), and the
		// background flow. Neither changes during the time step, so they
		// are calculated once, before the sub-steps.
		const double s_coeff {nanbu_s_coeff_cuda(q, T, Te, ne, slots_d.mass,
			mass_kg)};
		const double uX {Interpolate::interpolate_field_4d(bkg_d.uX, 
			stencil_4d, t, x, y, z)};
		const double uY {Interpolate::interpolate_field_4d(bkg_d.uY, 
			stencil_4d, t, x, y, z)};
		const double uZ {Interpolate::interpolate_field_4d(bkg_d.uZ, 
			stencil_4d, t, x, y, z)};

		// Sub-cycling. Nanbu's model collides the impurity with ONE
		// background particle per call. The friction force comes mostly
		// from slow background particles (s ~ 1 / g^3), and when the time
		// step is long for them (s >> 1) that one collision is already
		// fully random: it can't count for more, though in reality the
		// impurity would have met many slow particles in that time. The
		// force then comes out too weak. So when the step is too
		// collisional, it is split into n_sub collisions of dt / n_sub,
		// each with a newly sampled background particle. With n_sub = 1
		// this is the plain model. See the notes in collisions.h.
		const int n_sub {nanbu_num_substeps_cuda(
			nanbu_s0_cuda(s_coeff, dt, T, mass_kg),
			opts_d->imp_coll_subcycle_target,
			opts_d->imp_coll_subcycle_cap)};
		const double dt_sub {dt / n_sub};

		for (int i_sub {}; i_sub < n_sub; ++i_sub)
		{

		// The Nanbu model has three main variables in it:
		// s:    How collisional is this (sub-)step?
		// A(s): What is the shape of the scattering distribution? A(s) is
		//       the PDF.
		// chi:  What is the actual deflection angle this time? It is calculated
		//       via direct inversion from the CDF formed from our PDF (A(s)). 

		// Sample a background particle and calculate s (Eq. 19) for it,
		// with the impurity velocity as it is after the previous sub-steps.
		NanbuSResult s_res {nanbu_calc_s_cuda(s_coeff, dt_sub, T, mass_kg,
			uX, uY, uZ, vX, vY, vZ, rngs_d[i])};
		double s {s_res.s};
		double gX {s_res.gX};
		double gY {s_res.gY};
		double gZ {s_res.gZ};

		// Add to running sum of s values in each cell so we can do an average
		// later. Only consider for ions since they are the dominant collision
		// but no reason this can't be expanded for electrons as well. 
		// Generally leave this commented out unless you are investigating the 
		// collision model.
		/*
		if (!elec)
		{
			int Nx {stats_d.Nx};
			int Ny {stats_d.Ny};
			int Nz {stats_d.Nz};
			int idx {slots_d.tidx[i] * (Nx * Ny * Nz) + slots_d.xidx[i] 
				* (Ny * Nz) + slots_d.yidx[i] * Nz + slots_d.zidx[i]};
			double p_w_dt {slots_d.weight[i] * dt_sub};  // Change to float

			atomicAdd(&stats_d.s[idx], s * p_w_dt);
		}
		*/

		// Calcluate A (Eq. 13)
		double A {nanbu_calc_A_cuda(s)};

		// Calculate deflection angle, chi (Eq. 17)
		double chi {nanbu_calc_chi_cuda(s, A, rngs_d[i])};
		if (isnan(chi))
		{
			printf("Error! chi = nan\n");
			printf("  T   = %.17e\n", T);
			printf("  Te  = %.17e\n", Te);
			printf("  ne  = %.17e\n", ne);
			printf("  chi = %.17e\n", chi);
			printf("  s   = %.17e\n", s);
			printf("  A   = %.17e\n", A);
		}

		// Calculate post-collision velocity (Eq. 20a)
		Vec3d v_post {nanbu_post_coll_cuda(gX, gY, gZ, chi, slots_d.mass, 
			mass_kg, vX, vY, vZ, rngs_d[i])};

		// The next sub-step collides from this velocity
		vX = v_post.x;
		vY = v_post.y;
		vZ = v_post.z;

		}  // sub-steps

		// Update impurity velocity, which is v_n+1/2 since this is happening
		// after the Boris update.
		slots_d.vX[i] = vX;
		slots_d.vY[i] = vY;
		slots_d.vZ[i] = vZ;

	}  // nanbu_coll_kernel


	/**
	* Launches nanbu_coll_kernel to update the velocity of every particle on
	* the GPU due to collisions with one background species.
	*
	* Uses blocks of 256 threads, one thread per slot. The launch is
	* asynchronous, except in DEBUG builds, where it waits for the kernel to
	* finish and prints any error.
	*
	* @param slots_d Particle slots on the device.
	* @param bkg_d Background plasma on the device.
	* @param elec If true, collide with electrons. Otherwise collide with
	*        ions.
	* @param dt Time step (s).
	* @param imp_stats_d Impurity statistics on the device.
	* @param rngs_d Random number generators on the device, one per slot.
	* @param elec_mass_amu Electron mass (amu).
	* @param ion_mass_amu Main ion mass (amu).
	* @param opts_d Simulation options on the device.
	*/
	void collision_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d, const bool elec, 
		const double dt, ImpurityStats::StatisticsDevice& imp_stats_d,
		pcg32* rngs_d, const double elec_mass_amu, const double ion_mass_amu,
		const Options::OptionsDevice* opts_d)
	{
		
		// Block and grid size
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;

		nanbu_coll_kernel<<<gridSize, blockSize>>>(slots_d, bkg_d, elec, dt, 
			imp_stats_d, rngs_d, elec_mass_amu, ion_mass_amu, opts_d);

#ifdef DEBUG
		// Check for errors
		cudaError_t err {cudaDeviceSynchronize()};
		if (err != cudaSuccess)
			printf("nanbu_coll_kernel error: %s\n", cudaGetErrorString(err));
#endif

	}  // collision_gpu

}  // namespace Collisions
