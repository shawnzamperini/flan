cmake_minimum_required(VERSION 3.26)
project(flan VERSION 0.12.0)

# ------------------------------------------------------------
# Build types (Debug, Release, RelWithDebInfo, MinSizeRel)
# ------------------------------------------------------------
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the build type")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ------------------------------------------------------------
# Load dependency options (FLANSOFT, FLANROOT)
# ------------------------------------------------------------
set(INSTALL_DEPS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/install-deps")
include("${INSTALL_DEPS_PATH}/build-opts.cmake" REQUIRED)

message(STATUS "FLANSOFT set to ${FLANSOFT}")
message(STATUS "FLANROOT set to ${FLANROOT}")

# ------------------------------------------------------------
# Git commit hash
# ------------------------------------------------------------
execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "Commit: ${GIT_COMMIT_HASH}")

# ------------------------------------------------------------
# Version header
# ------------------------------------------------------------
message(STATUS "Flan version: ${flan_VERSION}")
configure_file(config.h.in config.h)
include_directories(${CMAKE_BINARY_DIR})

# ------------------------------------------------------------
# Language standards
# ------------------------------------------------------------
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# ------------------------------------------------------------
# Optional CUDA support
# ------------------------------------------------------------
option(ENABLE_CUDA "Enable CUDA GPU acceleration" OFF)

if(ENABLE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)

    # Needed for atomicAdd(double*, double). Requires CUDA Toolkit 12, 13
	# abandonded support for these? Should address...
    set(CMAKE_CUDA_ARCHITECTURES 60 70 80 86)

    message(STATUS "NVIDIA GPU support enabled")
else()
    message(STATUS "NVIDIA GPU support disabled")
endif()

# ------------------------------------------------------------
# MPI + OpenMP
# ------------------------------------------------------------
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

# ------------------------------------------------------------
# Source files
# ------------------------------------------------------------
file(GLOB CPU_SOURCES "src/*.cpp")

if(ENABLE_CUDA)
    file(GLOB CUDA_SOURCES "src/cuda/*.cu")
endif()

set(SOURCES ${CPU_SOURCES} ${CUDA_SOURCES})

# ------------------------------------------------------------
# Include directories
# ------------------------------------------------------------
include_directories(
    include
    ${FLANSOFT}/netcdf-cxx4/include
    ${FLANSOFT}/netcdf-c/include
)

if(ENABLE_CUDA)
    include_directories(include/cuda)
endif()

# ------------------------------------------------------------
# Link directories
# ------------------------------------------------------------
link_directories(
    ${FLANSOFT}/netcdf-cxx4/lib
    ${FLANSOFT}/netcdf-cxx4/lib64
    ${FLANSOFT}/netcdf-c/lib
)

# ------------------------------------------------------------
# Linux shared library flags
# ------------------------------------------------------------
if(UNIX AND NOT APPLE)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
endif()

# ------------------------------------------------------------
# Build library
# ------------------------------------------------------------
add_library(flan SHARED ${SOURCES})

# ------------------------------------------------------------
# CPU compile flags (per build type)
# ------------------------------------------------------------
target_compile_options(flan PRIVATE
    $<$<CONFIG:Release>:-O3 -Wall -fopenmp -march=native>
    $<$<CONFIG:RelWithDebInfo>:-O2 -g -Wall -fopenmp -march=native>
    $<$<CONFIG:Debug>:-O0 -g -Wall -fopenmp -march=native>
)

# ------------------------------------------------------------
# CUDA compile flags (per build type)
# ------------------------------------------------------------
if(ENABLE_CUDA)
    target_compile_definitions(flan PUBLIC USE_CUDA)
    target_link_libraries(flan PRIVATE CUDA::cudart)

    target_compile_options(flan PRIVATE
        $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CUDA>:-arch=sm_86>>
        $<$<CONFIG:RelWithDebInfo>:$<$<COMPILE_LANGUAGE:CUDA>:-arch=sm_86 -lineinfo>>
        $<$<CONFIG:Debug>:$<$<COMPILE_LANGUAGE:CUDA>:-G -O0>>
    )

    set_target_properties(flan PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()

# ------------------------------------------------------------
# Link libraries
# ------------------------------------------------------------
target_link_libraries(flan
    PRIVATE
        netcdf
        netcdf-cxx4
        OpenMP::OpenMP_CXX
)

# MPI must be public, since the input file (the executable) needs to call
# MPI_Init(). 
target_link_libraries(flan PUBLIC MPI::MPI_CXX)

# ------------------------------------------------------------
# RPATH for netCDF
# ------------------------------------------------------------
set_target_properties(flan PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${FLANSOFT}/netcdf-c/lib;${FLANSOFT}/netcdf-cxx4/lib;${FLANSOFT}/netcdf-cxx4/lib64"
)

# ------------------------------------------------------------
# Install library
# ------------------------------------------------------------
install(TARGETS flan
    EXPORT flan_targets
    LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
    PUBLIC_HEADER DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# ------------------------------------------------------------
# Export config
# ------------------------------------------------------------
install(EXPORT flan_targets
    FILE flanConfig.cmake
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/cmake/flan
)

# ------------------------------------------------------------
# Generate flan_setup.sh
# ------------------------------------------------------------
add_custom_command(
    OUTPUT script_output.txt
    COMMAND bash ${FLANROOT}/scripts/generate_flan_setup.sh ${FLANROOT} "${INSTALL_DEPS_PATH}/build-opts.sh"
    WORKING_DIRECTORY ${FLANROOT}/scripts
    COMMENT "Generating flan_setup.sh"
)

add_custom_command(
    OUTPUT flan_setup_copied
    DEPENDS script_output.txt
    COMMAND ${CMAKE_COMMAND} -E copy
            ${FLANROOT}/scripts/flan_setup.sh
            ${CMAKE_BINARY_DIR}/flan_setup.sh
    COMMENT "Copying flan_setup.sh to build directory"
)

add_custom_target(flan_setup ALL DEPENDS flan_setup_copied)

# ------------------------------------------------------------
# Enable unit tests
# ------------------------------------------------------------

# First get GoogleTest
include(FetchContent)
FetchContent_Declare(googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip)

# For MSVC: avoid overriding runtime library settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(INSTALL_GMOCK OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Enable tests
enable_testing()
add_subdirectory(tests)

