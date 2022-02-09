#!/bin/bash

module load 2019
module load GCC/7.3.0-2.30 OpenMPI/4.0.3-GCC-7.3.0-2.30 CMake/3.12.1-GCCcore-7.3.0

export MUSCLE3_HOME="$HOME/muscle3"

export OMPI_CXX="g++"
export OMPI_CC="gcc"
export MPI_COMPILER="mpicxx"
export MPI_COMPILER_FLAGS="-Wall -Wnon-virtual-dtor"
export MPI_LINKER_FLAGS=""

./build.sh "$@"
