#!/bin/bash

module load pre2019
module load OpenMPI/2.1.1-GCC-6.4.0-2.28
module load CMake/3.12.1-GCCcore-6.4.0

export MUSCLE3_HOME="$HOME/muscle3"

export OMPI_CXX="g++"
export OMPI_CC="gcc"
export MPI_COMPILER="mpicxx"
export MPI_COMPILER_FLAGS="-Wall -Wnon-virtual-dtor"
export MPI_LINKER_FLAGS=""

./build.sh "$@"
