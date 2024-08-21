#!/bin/bash

# Source the paths used for this installation
source build-opts.sh

# Install prefix
PREFIX=$FLANSOFT/postgkyl

# Delete old directory if it exists
rm -rf $PREFIX

# The python environment needed for Flan requires postgkyl. postgkyl is not
# available via conda or pip, only through github. Therefore we setup our
# python environment as so:
# 1. Clone postgkyl repository
git clone https://github.com/ammarhakim/postgkyl.git $PREFIX 
cd $PREFIX

# 2. Checkout commit from 8/20/24. We can update (or comment this line out) 
# for a more recent commit, but we're doing this just to keep things 
# constant.
git -c advice.detachedHead=false checkout b11210b

# 3. Create the flan conda environment using the postgkyl one as a starting
# point. Delete the flan environment if it already exists.
eval "$(conda shell.bash hook)"
conda env remove -n flan -y
conda env create -f environment.yml -n flan
conda activate flan

# 4. Install the environment with pip so we can find postgkyl without needing
# to modify PYTHONPATH and all that jank.
pip install -e .

# 5. Install additional dependencies needed by flan here
conda install conda-forge::netcdf4 -y
conda install ipython -y

# 6. Return back to the top flan directory and install the python scripts as
# the flan package so they can be found by the main code (and future python
# scripts as well).
cd $FLANROOT
pip install -e . 

# 7. Flan is coupled to postgkyl via the script python/read_gkyl.py. Therefore,
# Flan needs to know where to find this script. This is implemented via an
# environment variable called READ_GKYL_PY. 
conda env config vars set READ_GKYL_PY=$FLANROOT/python/read_gkyl.py
