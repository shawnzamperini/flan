#include "slots_device.h"


// Initializes RNGs on the GPU 
pcg32* init_rngs_cuda(Slots::SlotsDevice& slots_d, uint64_t base_seed);
