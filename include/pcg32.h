#pragma once
#include <cmath>
#include <cstdint>

#include "constants.h"

// If using nvcc define this so that this struct can be used in CPU and GPU code
#ifdef __CUDACC__
#define HD __host__ __device__
#else
#define HD
#endif

// I just had AI give me this, very straightforward RNG that is STL-free and
// thus can be used in CUDA
struct pcg32 
{
    uint64_t state;
    uint64_t inc;

	bool has_spare {false};
	double spare {};

	HD
    pcg32(uint64_t seed, uint64_t seq) {
        state = 0;
        inc = (seq << 1u) | 1u;
        next_uint();
        state += seed;
        next_uint();
    }

	HD
    uint32_t next_uint() {
        uint64_t oldstate = state;
        state = oldstate * 6364136223846793005ULL + inc;
        uint32_t xorshifted = static_cast<uint32_t>(((oldstate >> 18u) ^ oldstate) >> 27u);
        uint32_t rot = static_cast<uint32_t>(oldstate >> 59u);
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    }

	HD
    double next_double() {
        return (next_uint() >> 8) * (1.0 / 16777216.0);
    }

	// Box-Mueller transform for returning two normally distributed random
	// numbers centered on mean with standard deviation stddev. Since it
	// returns two numbers, we cache the second one to avoid throwing it away.
	HD
	double normal(double mean, double stddev)
	{
		// If spare random number was cached use it
		if (has_spare)
		{
			has_spare = false;
			return mean + stddev * spare;
		}

		// Box-Mueller transform for two normally distributed random numbers
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
