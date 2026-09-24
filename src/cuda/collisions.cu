#include <cuda_runtime.h>
#include <cstdio>
#include <math.h>

#include "background_device.h"
#include "constants.h"
#include "pcg32.h"
#include "slots_device.h"

#include "impurity_stats.cuh"
#include "interpolate.cuh"
#include "nanbu_s_a.cuh"
#include "utilities.cuh"


namespace Collisions
{
	// Simple struct to hold values returned from calc_s.
	struct NanbuSResult
	{
		double s;
		double gX;
		double gY;
		double gZ;
	};

	// Simple struct to hold three floats
	struct Vec3d
	{
		double x;
		double y;
		double z;
	};


	/**
	* Computes the electron Debye length.
	*
	* λ_D = sqrt(eps0 kT_e / (n_e e²))
	*
	* @param te Electron temperature [eV].
	* @param ne Electron density [m-3].
	*
	* @return Debye length [m].
	*/
	__device__ __forceinline__
	double calc_debye_length_cuda(const double te, const double ne)
	{
		return sqrt(Constants::eps0 * te * Constants::ev_to_j / 
			(ne * Constants::charge_e * Constants::charge_e)); 
	}


	/**
	 * Samples a velocity from a drifting Maxwellian distribution.
	 *
	 * The distribution has temperature T and bulk velocity
	 * (uX, uY, uZ).
	 *
	 * @param T Temperature [eV].
	 * @param uX Drift velocity x-component [m/s].
	 * @param uY Drift velocity y-component [m/s].
	 * @param uZ Drift velocity z-component [m/s].
	 * @param m Particle mass [kg].
	 * @param rng Random number generator.
	 * @param[out] bkg_vX Sampled x-velocity [m/s].
	 * @param[out] bkg_vY Sampled y-velocity [m/s].
	 * @param[out] bkg_vZ Sampled z-velocity [m/s].
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


	// Calculate s term in Nanbu collision model
	__device__ __forceinline__
	NanbuSResult nanbu_calc_s_cuda(const Background::BackgroundDevice& bkg_d,
		const double dt, const double T, const double mass_kg, const double Te,
		const double ne, const int tidx, const int xidx, const int yidx, 
		const int zidx, const double t, const double x, const double y,
		const double z, const double vX, const double vY, const double vZ,
		const int q, const double imp_mass, pcg32& rng, 
		Interpolate::InterpolationStencil stencil_4d)
	{

		// Interpolate to get background flow at particle location
		const double uX {Interpolate::interpolate_field_4d(bkg_d.uX, 
			stencil_4d, t, x, y, z)};
		const double uY {Interpolate::interpolate_field_4d(bkg_d.uY, 
			stencil_4d, t, x, y, z)};
		const double uZ {Interpolate::interpolate_field_4d(bkg_d.uZ, 
			stencil_4d, t, x, y, z)};

		// Sample background particle velocity
		double bkg_vX, bkg_vY, bkg_vZ;
		sample_bkg_velocity_cuda(T, uX, uY, uZ, mass_kg, rng, bkg_vX, bkg_vY, 
			bkg_vZ);

		// Instantaneous relative velocity
		double inst_gX = vX - bkg_vX;
		double inst_gY = vY - bkg_vY;
		double inst_gZ = vZ - bkg_vZ;

		// Mean relative velocity. Not actually part of the Nanbu algorithm, 
		// but if you want to match analytic fluid results this is often
		// a built in assumption in those derivations, so you'd swap inst_g
		// out with mean_g. inst_g is more kinetic/accurate though.
		const double mean_gX = vX - uX;
		const double mean_gY = vY - uY;
		const double mean_gZ = vZ - uZ;

		// Magnitude
		double mean_g = sqrt(mean_gX * mean_gX + mean_gY * mean_gY +
				 mean_gZ * mean_gZ);
		double inst_g = sqrt(inst_gX * inst_gX + inst_gY * inst_gY +
				 inst_gZ * inst_gZ);

		mean_g = fmax(mean_g, 1.0e-3);
		inst_g = fmax(inst_g, 1.0e-3);

		// Optional correction left disabled as in CPU version
		//inst_gX *= mean_g / inst_g;
		//inst_gY *= mean_g / inst_g;
		//inst_gZ *= mean_g / inst_g;

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

		// Calculate s
		const double square_term = q * Constants::charge_e * 
			Constants::charge_e / (Constants::eps0 * mu_ab);
		//const double s = ln_alpha / (4.0 * Constants::pi) * square_term *
		//	square_term * ne / (inst_g * inst_g * inst_g) * dt;

		// Needed to match assumptions made in the equation for the friction 
		// force
		const double s = ln_alpha / (4.0 * Constants::pi) * square_term *
			square_term * ne / (mean_g * mean_g * mean_g) * dt;

		// Return as struct
		//return {s, inst_gX, inst_gY, inst_gZ};
		return {s, mean_gX, mean_gY, mean_gZ};

	}  // nanbu_calc_s_cuda


	// Calculate A term in Nanbu collision model
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
			return Utilities::linear_interpolate_cuda(Nanbu::s, Nanbu::A, s);
		}
	}  // nanbu_calc_A


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

	
	// Calculate post-collision Cartesian velocity components from Nanbu model
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


	__global__
	void nanbu_coll_kernel(Slots::SlotsDevice slots_d, 
		const Background::BackgroundDevice bkg_d, bool elec, const double dt, 
		ImpurityStats::StatisticsDevice stats_d, pcg32* rngs_d, 
		const double elec_mass_amu, const double ion_mass_amu)
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
		int q {slots_d.q[i]};

		// Interpolate for ne, Te
		// Get interpolation stencil at this location for 4D arrays
		Interpolate::InterpolationStencil stencil_4d
			{Interpolate::build_stencil_4d(bkg_d, tidx, xidx, yidx, zidx)};

		// Interpolate to get background flow at particle location
		const double ne {Interpolate::interpolate_field_4d(bkg_d.ne, 
			stencil_4d, t, x, y, z)};
		const double Te {Interpolate::interpolate_field_4d(bkg_d.te, 
			stencil_4d, t, x, y, z)};

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
			T = fmax(0.1, T);
		}

		// The Nanbu model has three main variables in it:
		// s:    How collisional is this step?
		// A(s): What is the shape of the scattering distribution? A(s) is
		//       the PDF.
		// chi:  What is the actual deflection angle this time? It is calculated
		//       via direct inversion from the CDF formed from our PDF (A(s)). 

		// Calculate s (Eq. 19), making sure to pass in the full time step
		// velocities.
		NanbuSResult s_res {nanbu_calc_s_cuda(bkg_d, dt, T, mass_kg, Te, ne, 
			tidx, xidx, yidx, zidx, t, x, y, z, vX, vY, vZ, q, slots_d.mass, 
			rngs_d[i], stencil_4d)};
		double s {s_res.s};
		double gX {s_res.gX};
		double gY {s_res.gY};
		double gZ {s_res.gZ};

		// Add to running sum of s values in each cell so we can do an average
		// later. Only consider for ions since they are the dominant collision
		// but no reason this can't be expanded for electrons as well. 
		// Generally leave this commented out unless you are investigating the 
		// collision model.
		if (!elec)
		{
			int Nx {stats_d.Nx};
			int Ny {stats_d.Ny};
			int Nz {stats_d.Nz};
			int idx {slots_d.tidx[i] * (Nx * Ny * Nz) + slots_d.xidx[i] 
				* (Ny * Nz) + slots_d.yidx[i] * Nz + slots_d.zidx[i]};
			double p_w {slots_d.weight[i]};  // Change to float
			atomicAdd(&stats_d.s[idx], s * p_w);
		}

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

		// Update impurity velocity, which is v_n+1/2 since this is happening
		// after the Boris update.
		slots_d.vX[i] = v_post.x;
		slots_d.vY[i] = v_post.y;
		slots_d.vZ[i] = v_post.z;

	}  // nanbu_coll


	void collision_gpu(Slots::SlotsDevice& slots_d, 
		const Background::BackgroundDevice& bkg_d, const bool elec, 
		const double dt, ImpurityStats::StatisticsDevice& imp_stats_d,
		pcg32* rngs_d, const double elec_mass_amu, const double ion_mass_amu)
	{
		
		// Block and grid size
		int blockSize = 256;
		int gridSize  = (slots_d.N + blockSize - 1) / blockSize;

		nanbu_coll_kernel<<<gridSize, blockSize>>>(slots_d, bkg_d, elec, dt, 
			imp_stats_d, rngs_d, elec_mass_amu, ion_mass_amu);

#ifdef DEBUG
		// Check for errors
		cudaError_t err {cudaDeviceSynchronize()};
		if (err != cudaSuccess)
			printf("nanbu_coll_kernel error: %s\n", cudaGetErrorString(err));
#endif

	} // nanbu_coll_gpu

}  // namespace Collisions
