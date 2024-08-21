#!/bin/bash

# Source the paths used for this installation
source build-opts.sh

# Install prefix
PREFIX=$FLANSOFT/netcdf-cxx4

# Delete old directory if it exists
rm -rf netcdf-cxx4

# We use github here since it looks like the tarballs are a bit out of date.
git clone https://github.com/Unidata/netcdf-cxx4.git
cd netcdf-cxx4
#autoreconf -if
#CC=$CC CXX=$CXX ./configure --prefix=$PREFIX
#make
#make check
#make install

mkdir build
cd build
CC=$CC CXX=$CXX cmake -DCMAKE_PREFIX_PATH="$FLANSOFT/netcdf-c;$FLANSOFT/hdf5" -DCMAKE_INSTALL_PREFIX=$PREFIX ..
make
make test
make install
