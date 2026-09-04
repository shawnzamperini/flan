# Build the test executables: CPU tests
add_executable(test_get_nearest_index_cpu test_get_nearest_index_cpu.cpp)
add_executable(test_get_nearest_cell_index_cpu test_get_nearest_cell_index_cpu.cpp)
add_executable(test_fill_slots_cpu test_fill_slots_cpu.cpp)
add_executable(test_boundary_cpu test_boundary_cpu.cpp)
add_executable(test_boris_cpu test_boris_cpu.cpp)
add_executable(test_openadas_cpu test_openadas_cpu.cpp)

# Build the test executables: CUDA tests
add_executable(test_neighbor_index_cuda test_neighbor_index_cuda.cpp)
add_executable(test_all_dead_cuda test_all_dead_cuda.cpp)
add_executable(test_boris_cuda test_boris_cuda.cpp)
add_executable(test_fill_slots_cuda test_fill_slots_cuda.cpp)
add_executable(test_openadas_cuda test_openadas_cuda.cpp)
add_executable(test_trilinear_cuda test_trilinear_cuda.cu)

# These need to include the definitons in the executable since they call kernels
add_executable(test_boundary_cuda test_boundary_cuda.cu 
	${PROJECT_SOURCE_DIR}/src/cuda/boundary.cu
	${PROJECT_SOURCE_DIR}/src/cuda/device_constants.cu)

# Test executables that call specific CUDA functions (instead of via a wrapper)
set_target_properties(test_trilinear_cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(test_boundary_cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Enable the tests to inherit the RPATH of the flan library so that it can
# find the needed shared libraries (such as libnetcdf.so).
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Link your main library + GoogleTest. Tests that involve OpenMP need to link
# against those as well as the flan library.
target_link_libraries(test_get_nearest_index_cpu PRIVATE flan gtest_main)
target_link_libraries(test_get_nearest_cell_index_cpu PRIVATE flan gtest_main)
target_link_libraries(test_fill_slots_cpu PRIVATE flan gtest_main OpenMP::OpenMP_CXX)
target_link_libraries(test_boundary_cpu PRIVATE flan gtest_main OpenMP::OpenMP_CXX)
target_link_libraries(test_boris_cpu PRIVATE flan gtest_main)
target_link_libraries(test_openadas_cpu PRIVATE flan gtest_main OpenMP::OpenMP_CXX)

target_link_libraries(test_fill_slots_cuda PRIVATE flan gtest_main OpenMP::OpenMP_CXX)
target_link_libraries(test_openadas_cuda PRIVATE flan gtest_main)
target_link_libraries(test_all_dead_cuda PRIVATE flan gtest_main)
target_link_libraries(test_boris_cuda PRIVATE flan gtest_main)
target_link_libraries(test_neighbor_index_cuda PRIVATE flan gtest_main)
target_link_libraries(test_boundary_cuda PRIVATE flan gtest_main)

# When testing the CUDA kernel directly, no wrapper, we need to link against
# CUDA truntime.
target_link_libraries(test_trilinear_cuda PRIVATE flan gtest_main CUDA::cudart)
target_link_libraries(test_boundary_cuda PRIVATE flan gtest_main CUDA::cudart)

# Register with CTest
add_test(NAME GetNearestIndexCPUTest COMMAND test_get_nearest_index_cpu)
add_test(NAME GetNearestCellIndexCPUTest COMMAND test_get_nearest_cell_index_cpu)
add_test(NAME FillSlotsCPUTest COMMAND test_fill_slots_cpu)
add_test(NAME BoundaryCPUTest COMMAND test_boundary_cpu)
add_test(NAME BorisCPUTest COMMAND test_boris_cpu)
add_test(NAME OpenADASCPUTest COMMAND test_openadas_cpu)

add_test(NAME BorisCUDATest COMMAND test_boris_cuda)
add_test(NAME FillSlotsCUDATest COMMAND test_fill_slots_cuda)
add_test(NAME OpenADASCUDATest COMMAND test_openadas_cuda)
add_test(NAME TrilinearCUDATest COMMAND test_trilinear_cuda)
add_test(NAME AllDeadCUDATest COMMAND test_all_dead_cuda)
add_test(NAME NeighborIndexCUDATest COMMAND test_neighbor_index_cuda)
add_test(NAME BoundaryCUDATest COMMAND test_boundary_cuda)
