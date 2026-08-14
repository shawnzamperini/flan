# Build the test executable
add_executable(test_get_nearest_index_cpu test_get_nearest_index_cpu.cpp)
add_executable(test_get_nearest_cell_index_cpu test_get_nearest_cell_index_cpu.cpp)
add_executable(test_fill_slots_cpu test_fill_slots_cpu.cpp)
add_executable(test_fill_slots_gpu test_fill_slots_gpu.cpp)
add_executable(test_boundary_cpu test_boundary_cpu.cpp)
add_executable(test_neighbor_index_cuda test_neighbor_index_cuda.cpp)
add_executable(test_all_dead_cuda test_all_dead_cuda.cpp)
add_executable(test_trilinear_cuda test_trilinear_cuda.cu)

# CUDA function tests
set_target_properties(test_trilinear_cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Ensure CUDA support is included
#if(ENABLE_CUDA)
#    target_compile_definitions(test_fill_slots_gpu PRIVATE USE_CUDA)
#    target_link_libraries(test_fill_slots_gpu PRIVATE CUDA::cudart flan)
#endif()


# Enable the tests to inherit the RPATH of the flan library so that it can
# find the needed shared libraries (such as libnetcdf.so).
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Link your main library + GoogleTest
target_link_libraries(test_get_nearest_index_cpu PRIVATE flan gtest_main)
target_link_libraries(test_get_nearest_cell_index_cpu PRIVATE flan gtest_main)
target_link_libraries(test_fill_slots_cpu PRIVATE flan gtest_main)
target_link_libraries(test_fill_slots_gpu PRIVATE flan gtest_main)
target_link_libraries(test_boundary_cpu PRIVATE flan gtest_main)
target_link_libraries(test_neighbor_index_cuda PRIVATE flan gtest_main)
target_link_libraries(test_all_dead_cuda PRIVATE flan gtest_main)

# When testing the CUDA kernel directly, no wrapper, we need to link against
# CUDA truntime.
target_link_libraries(test_trilinear_cuda PRIVATE flan gtest_main CUDA::cudart)

#include_directories("/home/zamp/github/flan/include")
#target_link_libraries(test_get_nearest_index_cpu PRIVATE flan MPI::MPI_CXX)
#target_link_libraries(test_get_nearest_cell_index_cpu PRIVATE flan MPI::MPI_CXX)
#target_link_libraries(test_fill_slots_cpu PRIVATE flan MPI::MPI_CXX)
#target_link_libraries(test_fill_slots_gpu PRIVATE flan MPI::MPI_CXX)
#target_link_libraries(test_boundary_cpu PRIVATE flan MPI::MPI_CXX)


#link_directories(
#    ${FLANSOFT}/netcdf-cxx4/lib
#    ${FLANSOFT}/netcdf-cxx4/lib64
#    ${FLANSOFT}/netcdf-c/lib
#)

# Register with CTest
add_test(NAME GetNearestIndexCpuTest COMMAND test_get_nearest_index_cpu)
add_test(NAME GetNearestCellIndexCpuTest COMMAND test_get_nearest_cell_index_cpu)
add_test(NAME FillSlotsCpuTest COMMAND test_fill_slots_cpu)
add_test(NAME FillSlotsGpuTest COMMAND test_fill_slots_gpu)
add_test(NAME BoundaryCpuTest COMMAND test_boundary_cpu)
add_test(NAME NeighborIndexCUDATest COMMAND test_neighbor_index_cuda)
add_test(NAME AllDeadCUDATest COMMAND test_all_dead_cuda)
add_test(NAME TrilinearCUDATest COMMAND test_trilinear_cuda)

