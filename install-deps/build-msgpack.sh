#!/bin/bash

# Source that paths used for this installation
source build-opts.sh

# Install prefix
PREFIX=$FLANSOFT/msgpack-c

# Delete old
#rm -rf msgpack-c-6.0.2.tar.gz

# Download and unpackhttps://github.com/msgpack/msgpack-c/releases/download/c-6.0.2/msgpack-c-6.0.2.tar.gz
#curl -L https://github.com/msgpack/msgpack-c/releases/download/c-6.0.2/msgpack-c-6.0.2.tar.gz > msgpack-c-6.0.2.tar.gz
#tar xvzf msgpack-c-6.0.2.tar.gz

git clone https://github.com/msgpack/msgpack-c.git
cd msgpack-c

# Set to a May 28, 2024 commit for version control
git checkout 820ccf1f1d919b98a0bbe1ed9e9a004e922a80fc
#git checkout cpp_master
cmake -DMSGPACK_CXX20=ON -DCMAKE_INSTALL_PREFIX=$PREFIX .  
cmake --build . --target install

# Build with CMake
#cd msgpack-c-6.0.2
#CC=$CC CXX=$CXX cmake -DMSGPACK_CXX20=ON -DCMAKE_INSTALL_PREFIX=$PREFIX .
#cmake --build . --target install
