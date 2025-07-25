#!/bin/bash

# Contents of flan_setup.sh
cat <<EOF > flan_setup.sh
#!/bin/bash

# This file is generated automatically when building Flan. It setups a 
# directory for a simulation with a name passed in via the command line. The
# required CMakeLists.txt file is put in the directory along with a blank
# input file. After adding any desired inputs to the input file, compile it
# with Flan and run. Example usage:
#
#  $ ./flan_setup simulation_name
#  $ cd simulation_name
#  $ cmake .
#
# Edit the file simulation_name.cpp with any desired inputs, then run:
#
#  $ make
#  $ ./simulation_name
#
# If you edit simulation_name.cpp, only make needs to be re-run.

# Make directory of simulation name and enter it
mkdir \$1 && cd \$1

# Copy over a blank input file and rename it, only if it doesn't already exist
if [ -e "\$1.cpp" ]; then
	echo "Error! \$1.cpp already exists in this directory."
else
	cp $1/regression/blank_input_file.cpp \$1.cpp
fi

# Create a simple CMakeLists.txt file
cat <<EOF > CMakeLists.txt
cmake_minimum_required(VERSION 3.26)
project(\$1 LANGUAGES CXX)
find_package(flan REQUIRED HINTS "$1/lib/cmake")
add_executable(\$1 \$1.cpp)
include_directories("$1/include")
target_link_libraries(\$1 PRIVATE flan)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
EOF

# Need to put an EOF at the end of flan_setup.sh. Couldn't figure out how to
# escape the previous EOF. If you know how, modify it. 
echo "EOF" >> flan_setup.sh

cat <<EOF >> flan_setup.sh

# Run CMake, setting up the Makefile for the user
cmake .
EOF

# Let the user execute it
chmod +x flan_setup.sh
