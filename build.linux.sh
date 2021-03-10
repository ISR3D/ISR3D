#!/bin/bash
export MUSCLE3_HOME="$HOME/apps/muscle3"

export OMPI_CXX="g++"
export OMPI_CC="gcc"
export MPI_COMPILER="mpicxx"
export MPI_COMPILER_FLAGS="-Wall -Wnon-virtual-dtor"
export MPI_LINKER_FLAGS=""

./build.sh "$@"
