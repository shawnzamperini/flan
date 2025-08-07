# This version is likely a bit too strict, it is just the one I built this with
cmake_minimum_required(VERSION 3.26)
project(flan VERSION 0.5.2)

# Load options from when the dependencies were built. Contains FLANSOFT
# and FLANROOT
include(${CMAKE_CURRENT_SOURCE_DIR}/install-deps/build-opts.cmake REQUIRED)

# User must specify where the dependencies are stored via the FLANSOFT variable
#if(NOT DEFINED FLANSOFT)
#	set(FLANSOFT $ENV{HOME}/flansoft)
#endif()
message(STATUS "FLANSOFT set to ${FLANSOFT}")
message(STATUS "FLANROOT set to ${FLANROOT}")

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
	${FLANSOFT}/netcdf-cxx4/lib64
	${FLANSOFT}/netcdf-c/lib
)

if(UNIX AND NOT APPLE) # This applies to Linux platforms
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
endif()

# Add a shared library
add_library(flan SHARED ${SOURCES})

# Compilation flags
target_compile_options(flan PRIVATE -Wall -O3 -fopenmp)

# Specify the libraries to link against. I am not fully certain why we need
# the -fopenmp flag here in the linking stage, but it's needed
target_link_libraries(flan PRIVATE -fopenmp netcdf netcdf-cxx4)

# Set the rpath so we can find the netCDF shared library
set_target_properties(flan PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${FLANSOFT}/netcdf-c/lib;${FLANSOFT}/netcdf-cxx4/lib;${FLANSOFT}/netcdf-cxx4/lib64"
)

# Install the library itself
install(TARGETS flan
    EXPORT flan_targets
    LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
	PUBLIC_HEADER DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Export the library's target configuration
install(EXPORT flan_targets
    FILE flanConfig.cmake
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/cmake/flan
)

# Generate the flan_setup.sh script that helps a user setup a simulation. 
# script_output.txt here isn't actually anything, just a placeholder so we 
# can execute the generate script as part of the build process.
add_custom_command(
	OUTPUT script_output.txt
	COMMAND bash ${FLANROOT}/scripts/generate_flan_setup.sh ${FLANROOT} ${FLANSOFT}
	WORKING_DIRECTORY ${FLANROOT}/scripts
	COMMENT "Generating flan_setup.sh"
)
add_custom_target(flan_setup ALL DEPENDS script_output.txt)
