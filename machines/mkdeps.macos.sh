# If we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export FLANSOFT=$HOME/flansoft
export MACHINE_NAME='macos'
cd install-deps

# 1) Build OpenMPI.
# ./mkdeps.sh CC=clang CXX=clang++ --build-openmpi=yes --build-luajit=no  --prefix=$GKYLSOFT

# 2) Build gkylzero and adios. MacOS can use clang++, but I'm using g++ instead to be consistent.
#./mkdeps.sh CXX=g++ MPICXX=$FLANSOFT/openmpi-4.0.5/bin/mpicxx --prefix=$FLANSOFT --build-adios=yes --build-openmpi=no
./mkdeps.sh CC=clang CXX=clang++ --prefix=$FLANSOFT --build-adios=yes
