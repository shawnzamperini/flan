#pragma once
#include <cstdint>

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
};
