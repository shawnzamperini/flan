# Build the test executable
add_executable(test_get_nearest_index_cpu test_get_nearest_index_cpu.cpp)
add_executable(test_get_nearest_cell_index_cpu test_get_nearest_cell_index_cpu.cpp)

# Link your main library + GoogleTest
target_link_libraries(test_get_nearest_index_cpu PRIVATE flan gtest_main)
target_link_libraries(test_get_nearest_cell_index_cpu PRIVATE flan gtest_main)

# Register with CTest
add_test(NAME GetNearestIndexCpuTest COMMAND test_get_nearest_index_cpu)
add_test(NAME GetNearestCellIndexCpuTest COMMAND test_get_nearest_cell_index_cpu)

