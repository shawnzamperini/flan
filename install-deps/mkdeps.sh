#!/bin/bash

# Defaults
PREFIX=$HOME/flansoft
FLANROOT="NOT DEFINED"

# default build options
CC=gcc
CXX=g++
MPICC=mpicc
MPICXX=mpicxx

# By default, don't build anything. Will check later to see if things
# should be installed.
BUILD_ADIOS2=
BUILD_OPENMPI=
BUILD_BOOST=
BUILD_MSGPACK=
BUILD_ZLIB=
BUILD_HDF5=
BUILD_NETCDF_C=
BUILD_NETCDF_CXX=
INSTALL_NANOFLANN=
SETUP_CONDA_ENV=

# By default make a new build-opts.sh file.
NEW_BUILD_OPTS="yes"

# ----------------------------------------------------------------------------
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------------

# Help

show_help() {
cat <<EOF

./mkdeps.sh CC=cc CXX=cxx MPICC=mpicc MPICXX=mpicxx

Build Gkyl dependencies. By default, only builds libraries that Gkyl 
needs that haven't yet been built or can't be found.

CC 
CXX
MPICC
MPICXX                      C, C++, MPI C and MPI C++ compilers to use

-h
--help                      This help.
--prefix=DIR                Prefix where dependencies should be installed.
                            Default is $HOME/flansoft
--flanroot=DIR              Flanroot is the top-level directory of Flan. This
                            should be taken care of automatically and the
                            user should nto need to specify it.
--boost-inc-dir             Directory where boost is installed (for msgpack)

The following flags specify which libraries to build. By default, only
builds libraries that haven't yet been built or can't be found. 
If you build libraries that depend on MPI please specify the MPI C 
and C++ compilers to use.

--build-gkylzero            Should we build Gkylzero?
--build-adios2              Should we build ADIOS?
--build-openmpi             Should we build OpenMPI?
--build-zlib                SHould we build zlib? Needed for HDF5 and NetCDF-C
--build-hdf5                Should we build HDF5? Needed for NetCDF-C
--build-netcdf-c            Should we build netCDF?
--build-netcdf-cxx          Should we build the netCDF C++ interface?
--install-nanoflann         Should we install the (header-only) nanoflann library?
--setup-conda-env           Should we setup the conda environment?
--new-build-opts            Should we create a new build-opts.sh file?

The behavior of the flags for library xxx is as follows:
--build-xxx=no              Don't build xxx, even if it can't be found in PREFIXDIR
--build-xxx=yes             Build xxx, no matter what
[no flag specified]         Build xxx only if it can't be found in PREFIXDIR (default)
EOF
}

# Helper functions

find_program() {
   prog=`command -v "$1" 2>/dev/null`
   if [ -n "$prog" ]
   then
      dirname "$prog"
   fi
}

die() {
   echo "$*"
   echo
   echo "Dependency builds failed."
   echo
   exit 1
}

# ----------------------------------------------------------------------------
# MAIN PROGRAM
# ----------------------------------------------------------------------------

# Parse options

while [ -n "$1" ]
do
   value="`echo $1 | sed 's/[^=]*.\(.*\)/\1/'`"
   key="`echo $1 | sed 's/=.*//'`"
   if `echo "$value" | grep "~" >/dev/null 2>/dev/null`
   then
      echo
      echo '*WARNING*: the "~" sign is not expanded in flags.'
      echo 'If you mean the home directory, use $HOME instead.'
      echo
   fi
   case "$key" in
   -h)
      show_help
      exit 0
      ;;
   --help)
      show_help
      exit 0
      ;;
   CC)
      [ -n "$value" ] || die "Missing value in flag $key."
      CC="$value"
      ;;
   CXX)
      [ -n "$value" ] || die "Missing value in flag $key."
      CXX="$value"
      ;;
   MPICC)
      [ -n "$value" ] || die "Missing value in flag $key."
      MPICC="$value"
      ;;
   MPICXX)
      [ -n "$value" ] || die "Missing value in flag $key."
      MPICXX="$value"
      ;;   
   --prefix)
      [ -n "$value" ] || die "Missing value in flag $key."
      PREFIX="$value"
      ;;
   --flanroot)
      [ -n "$value" ] || die "Missing value in flag $key."
      FLANROOT="$value"
      ;;
   --build-openmpi)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_OPENMPI="$value"
      ;;
   --build-adios2)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_ADIOS2="$value"
      ;;
   --build-boost)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_BOOST="$value"
      ;;
   --build-msgpack)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_MSGPACK="$value"
      ;;
   --build-zlib)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_ZLIB="$value"
      ;;
   --build-hdf5)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_HDF5="$value"
      ;;
   --build-netcdf-c)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_NETCDF_C="$value"
      ;;
   --build-netcdf-cxx)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_NETCDF_CXX="$value"
      ;;
   --install-nanoflann)
      [ -n "$value" ] || die "Missing value in flag $key."
      INSTALL_NANOFLANN="$value"
      ;;
   --new-build-opts)
      [ -n "$value" ] || die "Missing value in flag $key."
      NEW_BUILD_OPTS="$value"
      ;;
   --setup-conda-env)
      [ -n "$value" ] || die "Missing value in flag $key."
      SETUP_CONDA_ENV="$value"
      ;;
   *)
      die "Error: Unknown flag: $1"
      ;;
   esac
   shift
done

# if mpicc doesn't work (because it doesn't exist or it's not in path), try to use installed openmpi version
if ! [ -x "$(command -v $MPICC)" ]
then
    MPICC=$PREFIX/openmpi/bin/mpicc
    MPICXX=$PREFIX/openmpi/bin/mpicxx
fi
# if mpicc still doesn't work, force to install openmpi
if ! [ -x "$(command -v $MPICC)" ] 
then
    BUILD_OPENMPI="yes"
fi

# Write out build options for scripts to use if specified
if [ "$NEW_BUILD_OPTS" = "yes" ]
then
cat <<EOF1 > build-opts.sh
# Generated automatically! Do not edit

# Installation directory
FLANSOFT=$PREFIX
FLANROOT=$FLANROOT
# Various compilers
CC=$CC
CXX=$CXX
MPICC=$MPICC
MPICXX=$MPICXX

EOF1
fi

build_adios2() {
    if [[ ! "$BUILD_ADIOS2" = "no" && ("$BUILD_ADIOS2" = "yes" || ! -f $PREFIX/adios/include/adios.h) ]]
    then    
	echo "Building ADIOS2"
	./build-adios2.sh
    fi
}

build_boost() {
   if [ "$BUILD_BOOST" = "yes" ]
   then
      echo "Building boost"
      ./build-boost.sh
   fi
}

build_msgpack() {
   if [ "$BUILD_MSGPACK" = "yes" ]
   then
      echo "Building msgpack"
      ./build-msgpack.sh
   fi
}

build_zlib() {
   if [ "$BUILD_ZLIB" = "yes" ]
   then
      echo "Building zlib"
      ./build-zlib.sh
   fi
}

build_hdf5() {
   if [ "$BUILD_HDF5" = "yes" ]
   then
      echo "Building HDF5"
      ./build-hdf5.sh
   fi
}

build_netcdf-c() {
   if [ "$BUILD_NETCDF_C" = "yes" ]
   then
      echo "Building netCDF-c"
      ./build-netcdf-c.sh
   fi
}

build_netcdf-cxx() {
   if [ "$BUILD_NETCDF_CXX" = "yes" ]
   then
      echo "Building netCDF-cxx"
      ./build-netcdf-cxx.sh
   fi
}

install_nanoflann() {
   if [ "$INSTALL_NANOFLANN" = "yes" ]
   then
      echo "Installing nanoflann"
      ./install-nanoflann.sh
   fi
}

setup_conda_env() {
   if [ "$SETUP_CONDA_ENV" = "yes" ]
   then
      echo "Setting up conda environment"
      ./setup-conda-env.sh
   fi
}

echo "Installations will be in $PREFIX"

# On the chopping block
#build_adios2
#build_boost
#build_msgpack

# Order matters for these four
build_zlib  # HDF5 (and thus netcdf-c) dependency
build_hdf5  # netcdf-c dependency
build_netcdf-c
build_netcdf-cxx

# Order not important for these
install-nanoflann
setup_conda_env
