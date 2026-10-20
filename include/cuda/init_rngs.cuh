#pragma once

/**
* @file  init_rngs.cuh
*
* @brief Header file for init_rngs.cu
*/
#include "pcg32.h"
#include "slots_device.h"


pcg32* init_rngs_cuda(Slots::SlotsDevice& slots_d, uint64_t base_seed);

void free_rngs(pcg32*& rngs_d, int device_id);
