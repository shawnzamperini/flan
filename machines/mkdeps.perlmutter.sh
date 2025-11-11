# If we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi

# We're now in the root directory of flan. Save to variable.
export FLANROOT=$PWD

# flansoft is where all the dependencies will live
export FLANSOFT=$HOME/flansoft
mkdir -p $FLANSOFT

export MACHINE_NAME='linux'
cd install-deps

# Build all the goodies. Specify Cray C/C++ compiler wrappers here.
./mkdeps.sh CC=cc CXX=CC --prefix=$FLANSOFT --flanroot=$FLANROOT --build-zlib=yes --build-hdf5=yes --build-netcdf-c=yes --build-netcdf-cxx=yes --setup-conda-env=yes
