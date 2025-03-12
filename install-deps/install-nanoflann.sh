#!/bin/bash

# Source the paths used for this installation
source build-opts.sh

# Install prefix
PREFIX=$FLANSOFT/nanoflann

# Delete old directory if it exists
rm -rf $PREFIX

# Clone nanoflann repository
git clone https://github.com/jlblancoc/nanoflann.git $PREFIX 
cd $PREFIX

# Checkout specific commit to keep things consistent. Can update if for
# whatever reason that is desired
git checkout 518b2b9

# Then just copy the heder file into the Flan include directory
cp include/nanoflann.hpp $FLANROOT/include
