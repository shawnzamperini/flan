#!/bin/bash

# Source the paths used for this installation
source build-opts.sh

# Install prefix
PREFIX=$FLANSOFT/hdf5

# Current directory, saving as variable since we need it later
INSTALL_DEPS_DIR=$(pwd)

# Delete old directory if it exists
rm -rf hdf5

# Download and unpack
curl -L https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.4/src/hdf5-1.14.4-3.tar.gz > hdf5-1.14.4-3.tar.gz
tar xvzf hdf5-1.14.4-3.tar.gz

# Make and install
cd hdf5-1.14.4-3
CC=$CC CXX=$CXX ./configure --with-zlib=$FLANSOFT/zlib --prefix=$PREFIX --enable-hl --enable-shared
make
make install

