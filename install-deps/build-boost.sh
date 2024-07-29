#!/bin/bash

# Source that paths used for this installation
source build-opts.sh

# Install prefix
PREFIX=$FLANSOFT/boost_1_85_0

# Delete old download
rm -rf boost_1_85_0.tar.gz boost_1_85_0
rm -rf $PREFIX/boost_1_85_0

# Download and unpack
curl -L https://archives.boost.io/release/1.85.0/source/boost_1_85_0.tar.gz > boost_1_85_0.tar.gz
tar xzf boost_1_85_0.tar.gz

# There's nothing that actually needs to be made since this is a header
# only library, so move it to the flansoft directory.
echo "Moving to $PREFIX"
mv boost_1_85_0 $PREFIX

# Save path to boost include directory so msgpack can find it
INC_LINE="Boost_INCLUDE_DIR=$PREFIX"
echo "Adding to build-opts: $INC_LINE"
grep -qxF "$INC_LINE" build-opts.sh || echo "$INC_LINE" >> build-opts.sh 
