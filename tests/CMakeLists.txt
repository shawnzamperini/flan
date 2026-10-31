# Enable the tests to inherit the RPATH of the flan library so that it can
# find the needed shared libraries (such as libnetcdf.so). This sets the
# default for targets created after it, so it has to come before them.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# ------------------------------------------------------------
# CPU tests. Built in every build.
# ------------------------------------------------------------

add_executable(test_get_nearest_index_cpu test_get_nearest_index_cpu.cpp)
add_executable(test_get_nearest_cell_index_cpu test_get_nearest_cell_index_cpu.cpp)
add_executable(test_fill_slots_cpu test_fill_slots_cpu.cpp)
add_executable(test_boundary_cpu test_boundary_cpu.cpp)
add_executable(test_boris_cpu test_boris_cpu.cpp)
add_executable(test_openadas_cpu test_openadas_cpu.cpp)
add_executable(test_zoned_vector_cpu test_zoned_vector_cpu.cpp)
add_executable(test_zones_cpu test_zones_cpu.cpp)

# Link your main library + GoogleTest. Tests that involve OpenMP need to link
# against those as well as the flan library.
target_link_libraries(test_get_nearest_index_cpu PRIVATE flan gtest_main)
target_link_libraries(test_get_nearest_cell_index_cpu PRIVATE flan gtest_main)
target_link_libraries(test_fill_slots_cpu PRIVATE flan gtest_main OpenMP::OpenMP_CXX)
target_link_libraries(test_boundary_cpu PRIVATE flan gtest_main OpenMP::OpenMP_CXX)
target_link_libraries(test_boris_cpu PRIVATE flan gtest_main)
target_link_libraries(test_openadas_cpu PRIVATE flan gtest_main OpenMP::OpenMP_CXX)
target_link_libraries(test_zoned_vector_cpu PRIVATE flan gtest_main)
target_link_libraries(test_zones_cpu PRIVATE flan gtest_main)

# Register with CTest
add_test(NAME GetNearestIndexCPUTest COMMAND test_get_nearest_index_cpu)
add_test(NAME GetNearestCellIndexCPUTest COMMAND test_get_nearest_cell_index_cpu)
add_test(NAME FillSlotsCPUTest COMMAND test_fill_slots_cpu)
add_test(NAME BoundaryCPUTest COMMAND test_boundary_cpu)
add_test(NAME BorisCPUTest COMMAND test_boris_cpu)
add_test(NAME OpenADASCPUTest COMMAND test_openadas_cpu)
add_test(NAME ZonedVectorCPUTest COMMAND test_zoned_vector_cpu)
add_test(NAME ZonesCPUTest COMMAND test_zones_cpu)

add_executable(test_cross_zones_cpu test_cross_zones_cpu.cpp)
target_link_libraries(test_cross_zones_cpu PRIVATE flan gtest_main
	OpenMP::OpenMP_CXX)
add_test(NAME CrossZonesCPUTest COMMAND test_cross_zones_cpu)

add_executable(test_zone_interp_cpu test_zone_interp_cpu.cpp)
target_link_libraries(test_zone_interp_cpu PRIVATE flan gtest_main
	OpenMP::OpenMP_CXX)
add_test(NAME ZoneInterpCPUTest COMMAND test_zone_interp_cpu)

add_executable(test_hdf5_utils_cpu test_hdf5_utils_cpu.cpp)
target_link_libraries(test_hdf5_utils_cpu PRIVATE flan gtest_main hdf5)
add_test(NAME Hdf5UtilsCPUTest COMMAND test_hdf5_utils_cpu)

add_executable(test_read_s3xe_cpu test_read_s3xe_cpu.cpp)
target_link_libraries(test_read_s3xe_cpu PRIVATE flan gtest_main hdf5)
add_test(NAME ReadS3XECPUTest COMMAND test_read_s3xe_cpu)

add_executable(test_s3x_bkg_cpu test_s3x_bkg_cpu.cpp)
target_link_libraries(test_s3x_bkg_cpu PRIVATE flan gtest_main hdf5
	OpenMP::OpenMP_CXX)
add_test(NAME S3XBkgCPUTest COMMAND test_s3x_bkg_cpu)

add_executable(test_stats_cpu test_stats_cpu.cpp)
target_link_libraries(test_stats_cpu PRIVATE flan gtest_main
	OpenMP::OpenMP_CXX)
add_test(NAME StatsCPUTest COMMAND test_stats_cpu)

add_executable(test_save_results_cpu test_save_results_cpu.cpp)
target_link_libraries(test_save_results_cpu PRIVATE flan gtest_main hdf5
	netcdf-cxx4 netcdf OpenMP::OpenMP_CXX)
add_test(NAME SaveResultsCPUTest COMMAND test_save_results_cpu)

# MPI test. It has its own main to initialize MPI, so it links gtest instead
# of gtest_main, and it's run on 3 ranks with the same mpiexec that CMake
# found for the library. The timeout stops ctest from waiting forever if a
# broken broadcast or reduction leaves the ranks hanging.
add_executable(test_mpi test_mpi.cpp)
target_link_libraries(test_mpi PRIVATE flan gtest hdf5)
add_test(NAME MPITest COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3
	${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi> ${MPIEXEC_POSTFLAGS})
set_tests_properties(MPITest PROPERTIES TIMEOUT 120)

# ------------------------------------------------------------
# CUDA tests. Only built with ENABLE_CUDA, since they need nvcc and the CUDA
# runtime (CUDA::cudart only exists once find_package(CUDAToolkit) has run).
# ------------------------------------------------------------

if(ENABLE_CUDA)

	add_executable(test_neighbor_index_cuda test_neighbor_index_cuda.cu)
	add_executable(test_all_dead_cuda test_all_dead_cuda.cpp)
	add_executable(test_boris_cuda test_boris_cuda.cpp)
	add_executable(test_fill_slots_cuda test_fill_slots_cuda.cpp)
	add_executable(test_openadas_cuda test_openadas_cuda.cpp)
	add_executable(test_trilinear_cuda test_trilinear_cuda.cu)
	add_executable(test_indices_cuda test_indices_cuda.cu)
	add_executable(test_step_cuda test_step_cuda.cpp)

	# These need to include the definitons in the executable since they call
	# kernels
	add_executable(test_boundary_cuda test_boundary_cuda.cu
		${PROJECT_SOURCE_DIR}/src/cuda/boundary.cu
		${PROJECT_SOURCE_DIR}/src/cuda/device_constants.cu)

	# So tricky. We have to compile with background.cu here so that we don't
	# use the flan library version of it, otherwise we end up with two
	# different versions of device_constants (d_x, etc.), which results in
	# undefined behavior. This occurs whenever to_device is used while needing
	# device_constants. Wish I didn't put myself in this situation, but unsure
	# how to fix it.
	add_executable(test_interpolate_cuda
		test_interpolate_cuda.cu
		${PROJECT_SOURCE_DIR}/src/cuda/device_constants.cu
		${PROJECT_SOURCE_DIR}/src/cuda/background.cu)

	# Test executables that call specific CUDA functions (instead of via a
	# wrapper)
	set_target_properties(test_trilinear_cuda PROPERTIES
		CUDA_SEPARABLE_COMPILATION ON)
	set_target_properties(test_boundary_cuda PROPERTIES
		CUDA_SEPARABLE_COMPILATION ON)
	set_target_properties(test_interpolate_cuda PROPERTIES
		CUDA_SEPARABLE_COMPILATION ON)

	target_link_libraries(test_fill_slots_cuda PRIVATE flan gtest_main OpenMP::OpenMP_CXX)
	target_link_libraries(test_openadas_cuda PRIVATE flan gtest_main CUDA::cudart)
	target_link_libraries(test_all_dead_cuda PRIVATE flan gtest_main)
	target_link_libraries(test_boris_cuda PRIVATE flan gtest_main)
	target_link_libraries(test_neighbor_index_cuda PRIVATE flan gtest_main)
	target_link_libraries(test_indices_cuda PRIVATE flan gtest_main)
	target_link_libraries(test_step_cuda PRIVATE flan gtest_main CUDA::cudart)
	target_link_libraries(test_boundary_cuda PRIVATE flan gtest_main)
	target_link_libraries(test_interpolate_cuda PRIVATE flan gtest_main)

	# When testing the CUDA kernel directly, no wrapper, we need to link
	# against CUDA runtime.
	target_link_libraries(test_trilinear_cuda PRIVATE flan gtest_main CUDA::cudart)
	target_link_libraries(test_boundary_cuda PRIVATE flan gtest_main CUDA::cudart)
	#target_link_libraries(test_interpolate_cuda PRIVATE flan gtest_main CUDA::cudart)

	# Register with CTest
	add_test(NAME BorisCUDATest COMMAND test_boris_cuda)
	add_test(NAME FillSlotsCUDATest COMMAND test_fill_slots_cuda)
	add_test(NAME OpenADASCUDATest COMMAND test_openadas_cuda)
	add_test(NAME TrilinearCUDATest COMMAND test_trilinear_cuda)
	add_test(NAME AllDeadCUDATest COMMAND test_all_dead_cuda)
	add_test(NAME NeighborIndexCUDATest COMMAND test_neighbor_index_cuda)
	add_test(NAME IndicesCUDATest COMMAND test_indices_cuda)
	add_test(NAME StepCUDATest COMMAND test_step_cuda)
	add_test(NAME BoundaryCUDATest COMMAND test_boundary_cuda)
	add_test(NAME InterpolateCUDATest COMMAND test_interpolate_cuda)

endif()
