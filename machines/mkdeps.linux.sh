# If we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi

# flansoft is where all the dependencies will live
export FLANSOFT=$HOME/flansoft
mkdir -p $FLANSOFT

export MACHINE_NAME='linux'
cd install-deps

# 1) Build OpenMPI.
# ./mkdeps.sh CC=clang CXX=clang++ --build-openmpi=yes --build-luajit=no  --prefix=$GKYLSOFT

# 2) Build gkylzero and adios. MacOS can use clang++, but I'm using g++ instead to be consistent.
#./mkdeps.sh CXX=g++ MPICXX=$FLANSOFT/openmpi-4.0.5/bin/mpicxx --prefix=$FLANSOFT --build-adios=yes --build-openmpi=no

# 1) Build boost first since it's needed for msgpack
#./mkdeps.sh CC=gcc CXX=g++ --prefix=$FLANSOFT --build-adios2=no --build-boost=yes --build-msgpack=no

# 2) Then msgpack now that we have boost and can tell it where to look
#./mkdeps.sh CC=gcc CXX=g++ --prefix=$FLANSOFT --build-adios2=no --build-boost=no --build-msgpack=yes --new-build-opts=no

# Build HDF5
#./mkdeps.sh CC=gcc CXX=g++ --prefix=$FLANSOFT --build-netcdf-c=no --build-adios2=no --build-msgpack=no --build-boost=no --build-hdf5=yes

# Build netCDF-c
./mkdeps.sh CC=gcc CXX=g++ --prefix=$FLANSOFT --build-adios2=no --build-msgpack=no --build-boost=no --build-zlib=yes --build-hdf5=yes --build-netcdf-c=yes
