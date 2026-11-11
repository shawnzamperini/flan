#!/bin/bash

# Source the paths used for this installation
source build-opts.sh

# Install prefix
PREFIX=$FLANSOFT/postgkyl

# Delete old directory if it exists
rm -rf $PREFIX

# The python environment needed for Flan requires postgkyl. postgkyl is not
# available via mamba or pip, only through github. Therefore we setup our
# python environment as so:
# 1. Clone postgkyl repository
git clone https://github.com/ammarhakim/postgkyl.git $PREFIX 
cd $PREFIX

# 2. Checkout commit from 8/20/24. We can update (or comment this line out) 
# for a more recent commit, but we're doing this just to keep things 
# constant.
git -c advice.detachedHead=false checkout b11210b

# 3. Create the flan mamba environment using the postgkyl one as a starting
# point. Delete the flan environment if it already exists. Must first
# initialize mamba and activate it since we're in a sub-shell.
eval "$(mamba shell hook --shell bash)"
mamba activate
mamba env remove -n flan -y
mamba env create -f environment.yml -n flan
mamba activate flan

# 4. Install the environment with pip so we can find postgkyl without needing
# to modify PYTHONPATH and all that jank.
pip install -e .

# 5. Install additional dependencies needed by flan here
mamba install -y conda-forge::netcdf4 ipython scipy tqdm \
	"conda-forge::rust>=1.85" conda-forge::h5py

# 6. Install RustBCA, which computes the sputtering yield tables used by
# Flan (python/rustbca_yields.py). Like postgkyl, it isn't on conda or pip,
# only github. It's written in Rust, and the Rust compiler (>= 1.85)
# installed into the environment in step 5 is used, so no system Rust is
# needed. The Python bindings are compiled on install, which takes a few
# minutes.
RUSTBCA_PREFIX=$FLANSOFT/RustBCA
rm -rf $RUSTBCA_PREFIX
git clone https://github.com/lcpp-org/RustBCA.git $RUSTBCA_PREFIX
cd $RUSTBCA_PREFIX

# Checkout commit from 10/19/26 (v3.0.0), to keep things constant.
git -c advice.detachedHead=false checkout 51bee889956147e3287491d0c68957967235a33c
echo "Building RustBCA, this takes a few minutes..."
PATH=$CONDA_PREFIX/bin:$PATH pip install .

# 7. Return back to the top flan directory and install the python scripts as
# the flan package so they can be found by the main code (and future python
# scripts as well).
cd $FLANROOT
pip install -e . 

# 8. Flan is coupled to postgkyl via the script python/read_gkyl.py. Therefore,
# Flan needs to know where to find this script. This is implemented via an
# environment variable called READ_GKYL_PY. 
conda env config vars set READ_GKYL_PY=$FLANROOT/python/read_gkyl.py

# 9. Again, but for the script used to calculate the electric field.
conda env config vars set CALC_ELEC_FIELD_PY=$FLANROOT/python/calc_elec_field.py

# 10. Check that RustBCA built. Only sputtering runs need it, so this is
# only a warning.
if ! python -c "import libRustBCA" 2>/dev/null; then
	echo "Warning! RustBCA (libRustBCA) didn't install. It's only needed to"
	echo "make sputtering yield tables. See the installation docs."
fi
