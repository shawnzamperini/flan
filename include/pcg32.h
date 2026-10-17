#pragma once

/**
* @file  pcg32.h
*
* @brief PCG32 random number generator that can be used in CPU and GPU code.
*
* PCG32 (the XSH RR variant of M. E. O'Neill's PCG family) has 64 bits of
* state and returns 32-bit numbers. It's small, fast and doesn't use the
* standard library, so the same generator works on both the host and the
* device. On the CPU there is one generator per OpenMP thread, and on the
* GPU one per slot (see init_rngs.cu). Generators with the same seed but
* different seq values (streams) give independent sequences.
*/
#include <cmath>
#include <cstdint>

#include "constants.h"

// If using nvcc define this so that this struct can be used in CPU and GPU
// code
#ifdef __CUDACC__
#define HD __host__ __device__
#else
#define HD
#endif

// I just had AI give me this, very straightforward RNG that is STL-free and
// thus can be used in CUDA
/**
* PCG32 random number generator, usable on the host and the device.
*/
struct pcg32 
{
    uint64_t state;	///< Internal state, advanced on every call
    uint64_t inc;	///< Stream increment, always odd. Set from seq.

	bool has_spare {false};	///< True if spare holds an unused normal value
	double spare {};		///< Second value from the last Box-Muller
							///< transform, returned by the next normal call

	/**
	* Creates a generator from a seed and a stream number.
	*
	* @param seed Starting seed.
	* @param seq Stream number. Generators with different seq values give
	*        independent sequences, even with the same seed.
	*/
	HD
    pcg32(uint64_t seed, uint64_t seq) {
        state = 0;
        inc = (seq << 1u) | 1u;
        next_uint();
        state += seed;
        next_uint();
    }

	/**
	* Returns the next random 32-bit unsigned integer.
	*
	* @return Uniformly distributed value in [0, 2^32 - 1].
	*/
	HD
    uint32_t next_uint() {
        uint64_t oldstate = state;
        state = oldstate * 6364136223846793005ULL + inc;
        uint32_t xorshifted = static_cast<uint32_t>(((oldstate >> 18u) ^ oldstate) >> 27u);
        uint32_t rot = static_cast<uint32_t>(oldstate >> 59u);
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    }

	/**
	* Returns the next random double in [0, 1).
	*
	* Uses two calls to next_uint to fill all 53 bits of a double's
	* mantissa, so the result is a multiple of 2^-53 (about 1e-16). With
	* fewer bits, comparisons like r < p would be wrong for small p, since
	* every p smaller than the spacing would act the same.
	*
	* @return Uniformly distributed value in [0, 1).
	*/
	HD
    double next_double() {
        // Top 27 bits of one number and top 26 bits of the next, combined
        // into a 53-bit integer and scaled by 2^-53
        const uint64_t hi = next_uint() >> 5;
        const uint64_t lo = next_uint() >> 6;
        return (hi * 67108864.0 + lo) * (1.0 / 9007199254740992.0);
    }

	/**
	* Returns a normally distributed random number.
	*
	* Uses the Box-Muller transform, which makes two independent normal
	* values from two uniform ones. The second is cached in spare and
	* returned by the next call, so only every other call draws new
	* uniform numbers. The cached value is in standard units, so the next
	* call can use a different mean and stddev.
	*
	* @param mean Mean of the distribution.
	* @param stddev Standard deviation of the distribution.
	* @return Random value from the normal distribution.
	*/
	HD
	double normal(double mean, double stddev)
	{
		// If spare random number was cached use it
		if (has_spare)
		{
			has_spare = false;
			return mean + stddev * spare;
		}

		// Box-Muller transform for two normally distributed random numbers.
		// u1 is kept above 0 so log(u1) is finite.
		double u1 = next_double();
		if (u1 < 1e-12) u1 = 1e-12;
		
		double u2 = next_double();

		double r = sqrt(-2.0 * log(u1));
		double theta = 2.0 * Constants::pi * u2;

		// Cache the second number
		spare = r * sin(theta);
		has_spare = true;

		return mean + stddev * (r * cos(theta));
	}
};
