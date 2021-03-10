#!/bin/bash

module load gcc/6.2.0 openmpi/3.1.4_gcc620

export MUSCLE3_HOME="$HOME/muscle3"

export OMPI_CXX="g++"
export OMPI_CC="gcc"
export MPI_COMPILER="mpicxx"
export MPI_COMPILER_FLAGS="-Wall -Wnon-virtual-dtor"
export MPI_LINKER_FLAGS=""

./build.sh "$@"
