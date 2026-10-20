#pragma once

/**
* @file  errors.h
*
* @brief Header file for errors.cpp
*/
#include <string>

#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif


namespace Errors
{
	[[noreturn]] void fatal_error(const std::string& msg);

#ifdef USE_CUDA
	void check_cuda(cudaError_t err, const char* context);
#endif
}
