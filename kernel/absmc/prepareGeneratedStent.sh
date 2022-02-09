#!/bin/bash

# This shell script runs the 1st three stages of ISR simulation: vessel generation, equilibration and stent deployment
# Argument is the name of .cfg file (usually in /cxa directory)

CONFIGNAME=$1
mkdir -p ../../cxa/data.equilibration
mkdir -p ../../cxa/data.deployment
./build/generateSmoothCylinder $CONFIGNAME
./build/equilibration3D $CONFIGNAME
./build/deployment3D-stent $CONFIGNAME
