# This version is likely a bit too strict, it is just the one I built this with
cmake_minimum_required(VERSION 3.31)
project(flan VERSION 0.1.0)

# User must specify where the dependencies are stored via the FLANSOFT variable
if(NOT DEFINED FLANSOFT)
	set(FLANSOFT $ENV{HOME}/flansoft)
endif()
message(STATUS "FLANSOFT set to ${FLANSOFT}")

# Print out version and make configuration header file so version is accessible
# within flan
message(STATUS "Flan version: ${flan_VERSION}")
configure_file(config.h.in config.h)

# Include directory where config.h is generated
include_directories(${CMAKE_BINARY_DIR})

# Enforce a relatively modern standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include all files in src
file(GLOB SOURCES "src/*.cpp")

# Add include directories for flan and the netcdf dependencies
include_directories(
    include
    ${FLANSOFT}/netcdf-cxx4/include
    ${FLANSOFT}/netcdf-c/include
)

# Link to the directories with each library
link_directories(
	${FLANSOFT}/netcdf-cxx4/lib
	${FLANSOFT}/netcdf-c/lib
)

# Add executable
add_executable(flan ${SOURCES})

# Compilation flags
target_compile_options(flan PRIVATE -Wall -O3 -fopenmp)

# Specify the libraries to link against. I am not fully certain why we need
# the -fopenmp flag here in the linking stage, but it's needed
target_link_libraries(flan PRIVATE -fopenmp netcdf netcdf-cxx4)

# Set the rpath so we can find the netCDF shared library
set_target_properties(flan PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH ${FLANSOFT}/netcdf-cxx4/lib
)
