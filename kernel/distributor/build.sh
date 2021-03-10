#!/bin/bash

# Shell script to build the kernel,
# basically a wrapper for cmake and make.


PROJECT_NAME="distributor"
BUILD_DIR="./build"
# run ccmake to configure a build of the stand-alone applications only

if [ $1 == "compile" ]; then
    if [[ ! -d $BUILD_DIR ]]; then mkdir $BUILD_DIR; fi
    cd $BUILD_DIR
    cmake ../src/ -DMUSCLE3_HOME=$MUSCLE3_HOME -DLOCAL_LIB=$ISR_LIB_DIR
    make -j8
fi

if [ $1 == "clean" ]; then
    echo "Cleaning $PROJECT_NAME"
    rm -rf $BUILD_DIR
fi

if [ $1 == "install" ]; then
    echo "Installing $PROJECT_NAME"
    cp $BUILD_DIR/distributor $ISR_BUILD_DIR
fi

