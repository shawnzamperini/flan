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

# Make and install
# --disable-byterange is if you do not have libcurl (having just curl is 
# not enough). It's easy to install this on a machine, but opting out
# just to avoid an extra step
CPPFLAGS="-I${FLANSOFT}/hdf5/include -I${FLANSOFT}/zlib/include" LDFLAGS="-L${FLANSOFT}/hdf5/lib -L${FLANSOFT}/zlib/lib" ./configure --prefix=$PREFIX --disable-byterange --with-hdf5=${FLANSOFT}/hdf5 --with-zlib=${FLANSOFT}/zlib --disable-dap
make
make check
make install

#echo "Deleting old build directory..."
#rm -rf build
#mkdir build
#cd build
#CC=$CC CXX=$CXX cmake -DCMAKE_PREFIX_PATH=$PREFIX -DENABLE_NETCDF_4=ON -DENABLE_DAP=OFF -DENABLE_HDF5=ON -DZLIB_INCLUDE_DIR=$FLANSOFT/zlib/include -DZLIB_LIBRARY=$FLANSOFT/zlib/lib -DHDF5_ROOT=$FLANSOFT/hdf5 ..

  ##
  # Accommodate developers who have hdf5 libraries and
  # headers on their system, but do not have a the hdf
  # .cmake files.  If this is the case, they should
  # specify HDF5_HL_LIBRARY, HDF5_LIBRARY, HDF5_INCLUDE_DIR manually.

#CC=$CC CXX=$CXX cmake -DCMAKE_PREFIX_PATH=$PREFIX -DENABLE_NETCDF_4=ON -DENABLE_DAP=OFF -DENABLE_HDF5=ON -DHDF5_ROOT=$FLANSOFT/hdf5 -DHDF5_HL_LIBRARY=$FLANSOFT/hdf5/lib/libhdf5_hl.so -DHDF5_LIBRARY=$FLANSOFT/hdf5/lib/libhdf5.so -DHDF5_INCLUDE_DIR=$FLANSOFT/hdf5/include ..
#make
#make test
#make install
#cmake --build . --target test
#cmake --build . --target install
