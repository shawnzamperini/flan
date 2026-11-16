# Builds the regression cases (reg_*.cpp) against an installed Flan.
#
#   cmake -S <flan>/regression -B <build dir> [-DFLAN_ROOT=<flan>]
#   cmake --build <build dir>
#
# FLAN_ROOT is the top directory of Flan, where "make install" put
# lib/cmake. It defaults to the directory above this one. If CMake picks up
# a different MPI than Flan was built with, pass the same compilers Flan
# was configured with, e.g. -DMPI_CXX_COMPILER=<path to mpicxx>.
cmake_minimum_required(VERSION 3.26)
project(flan_regression LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(FLAN_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." CACHE PATH
	"Top directory of Flan")
find_package(flan REQUIRED HINTS "${FLAN_ROOT}/lib/cmake")
find_package(MPI REQUIRED)

# One program per case
foreach(case reg_friction reg_thermal_force reg_drifts)
	if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${case}.cpp")
		add_executable(${case} ${case}.cpp)
		target_include_directories(${case} PRIVATE "${FLAN_ROOT}/include")
		target_link_libraries(${case} PRIVATE flan MPI::MPI_CXX)
	endif()
endforeach()
