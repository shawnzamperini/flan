# Our main controling executable
add_executable(flan flan_main.cpp read_gkyl.cpp)

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${HOME}/flansoft/adios2/lib)
find_library(LIB_ADIOS2 NAMES adios2_cxx11 PATHS ${HOME}/flansoft/adios2/lib)
target_link_libraries(flan ${LIB_ADIOS2}) 

# These lines are so we can include ADIOS2 (it depends on MPI)
#find_package(MPI REQUIRED)
#find_package(ADIOS2 REQUIRED)
#target_link_libraries(flan PRIVATE adios2::cxx11_mpi MPI::MPI_C)
#target_link_libraries(flan PRIVATE adios2::cxx11)

# Include path to the header files. 
target_include_directories(flan PRIVATE ${FLAN_ROOT}/include ${ADIOS_DIR}/include)

# Need to provide paths to synamic libraries that are loaded during runtime. 

