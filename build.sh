#!/bin/bash

#common build script. Instead of running this directly, use machine-specific scripts to set up compilers and MUSCLE3 location.

# get location of this script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

#directory names for all submodels in /kernel
moduleDirs="absmc voxelizer collector distributor flow3d"

export ISR_HOME="$DIR"
export ISR_LIB_DIR="$ISR_HOME/lib/"
export ISR_BUILD_DIR="$ISR_HOME/build/"
export PALABOS_DIR="$ISR_HOME/lib/palabos-v2.2.0"

echo "ISR3D located at" $ISR_HOME

if [ "$#" -lt 1 ] || [ $1 == "install" ]; then
    mkdir -p $ISR_BUILD_DIR
    for val in $moduleDirs; do
        cd "$ISR_HOME/kernel/$val"
        echo
        ./build.sh compile
        ./build.sh install
    done
    cd "$ISR_HOME"
fi


if [ "$#" -gt 0 ] && [ $1 == "clean" ]; then
    for val in $moduleDirs; do
        cd "$ISR_HOME/kernel/$val"
        ./build.sh clean
    done
    cd "$ISR_HOME"
fi
