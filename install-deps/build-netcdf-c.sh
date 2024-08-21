#!/bin/bash

# Source the paths used for this installation
source build-opts.sh

# Install prefix
PREFIX=$FLANSOFT/netcdf-c

# Delete old directory if it exists
rm -rf netcdf-c

# Download and unpack
curl -L https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz > netcdf-c-4.9.2.tar.gz
tar xvzf netcdf-c-4.9.2.tar.gz
cd netcdf-c-4.9.2

# Make and install. Need to include and link to HDF5 and zlib with CPPFLAGS,
# LDFLAGS and --with-hdf5 and --with-zblib. 
# --disable-byterange is if you do not have libcurl (having just curl is 
# not enough). It's easy to install on most machines, but opting out
# just to avoid an extra step.
# --disable-dap is some network thing that would add an unneeded dependency.
CC=$CC CXX=$CXX CPPFLAGS="-I${FLANSOFT}/hdf5/include -I${FLANSOFT}/zlib/include" LDFLAGS="-L${FLANSOFT}/hdf5/lib -L${FLANSOFT}/zlib/lib" ./configure --prefix=$PREFIX --disable-byterange --with-hdf5=${FLANSOFT}/hdf5 --with-zlib=${FLANSOFT}/zlib --disable-dap
make
#make check  # Optional
make install
