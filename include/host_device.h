#pragma once
/**
* @file host_device.h
*
* @brief The FLAN_HD macro, for functions shared by the CPU and GPU code.
*
* A function marked FLAN_HD is compiled for both the host and the device
* when nvcc compiles the file (a .cu file), and as an ordinary host
* function otherwise (a .cpp file, in a CPU or a GPU build). Writing a
* per-particle routine once this way keeps the CPU and GPU versions from
* drifting apart, since there is only one.
*
* Shared functions should be inline (or templates) in headers, take plain
* values, pointers and sizes rather than std::vector or the host classes,
* and not throw, allocate or print. See "CPU and GPU code" in the
* developer docs (code-structure.md).
*
* The test is on __CUDACC__ (set by nvcc), not USE_CUDA: in a GPU build,
* the .cpp files are still compiled by the host compiler, which doesn't
* know __host__ or __device__.
*/

#ifdef __CUDACC__
#define FLAN_HD __host__ __device__
#else
#define FLAN_HD
#endif
