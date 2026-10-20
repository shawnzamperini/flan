/**
* @file  errors.cpp
*
* @brief Error handling for errors that must stop the whole run.
*
* These are host-only. Device code can't end the program, so a kernel that
* finds an impossible state should call __trap(). That makes the next CUDA
* call on the host fail, and check_cuda turns the failure into a
* fatal_error. Problems a run can continue past should instead be counted
* on the device and reported by the host afterward (see the warnings in
* openadas.cu).
*/
#include <cstdlib>
#include <iostream>
#include <string>

#include "errors.h"
#include "mpi.h"

#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif


namespace Errors
{

	/**
	* Prints an error message and stops the whole run.
	*
	* Use this for errors that would give wrong results if Flan kept going.
	* The message is printed to stderr with the MPI rank, then MPI_Abort
	* ends every rank. std::exit would only end this rank, and the others
	* would hang at the next collective. If MPI isn't initialized (e.g., in
	* the unit tests), std::exit is used instead, since MPI_Abort can't be
	* called before MPI_Init.
	*
	* @param msg Description of the error. "Error! " is added to the front
	*        and a newline to the end.
	*/
	[[noreturn]] void fatal_error(const std::string& msg)
	{
		int mpi_initialized {};
		MPI_Initialized(&mpi_initialized);

		if (mpi_initialized)
		{
			int rank {};
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);

			// Write the message in one call so it doesn't get mixed in
			// with messages from other ranks
			std::cerr << ("Error (rank " + std::to_string(rank) + ")! "
				+ msg + '\n');
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}
		else std::cerr << ("Error! " + msg + '\n');

		// MPI_Abort isn't marked noreturn, so exit here either way
		std::exit(EXIT_FAILURE);
	}

#ifdef USE_CUDA

	/**
	* Stops the run if a CUDA call failed.
	*
	* Wrap CUDA runtime calls, or the cudaDeviceSynchronize after a kernel
	* launch, with this. For example:
	*
	*     check_cuda(cudaDeviceSynchronize(), "fill_slots_kernel");
	*
	* @param err Return value of the CUDA call.
	* @param context What was being done, for the error message (e.g., the
	*        kernel or function name).
	*/
	void check_cuda(cudaError_t err, const char* context)
	{
		if (err != cudaSuccess)
		{
			fatal_error(std::string(context) + ": "
				+ cudaGetErrorString(err));
		}
	}

#endif

}
