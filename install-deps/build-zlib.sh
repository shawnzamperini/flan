#!/bin/bash

# Source the paths used for this installation
source build-opts.sh

# Install prefix
PREFIX=$FLANSOFT/zlib

# Delete old directory if it exists
rm -rf zlib

# Download and unpack
curl -L https://www.zlib.net/zlib-1.3.1.tar.gz > zlib-1.3.1.tar.gz
tar xvzf zlib-1.3.1.tar.gz

# Make and install
cd zlib-1.3.1
CC=$CC CXX=$CXX ./configure --prefix=$PREFIX
make
make install

