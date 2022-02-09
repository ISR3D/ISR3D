#!/bin/bash

# This shell script runs the 1st three stages of ISR simulation: vessel generation, equilibration and stent deployment
# Argument is the name of .cfg file (usually in /cxa directory)

CONFIGNAME=$1
mkdir ./data.equilibration
mkdir ./data.deployment
/projects/0/einf462/ISR3D_UQ/ISR3D/kernel/absmc/build/generateSmoothCylinder $CONFIGNAME
/projects/0/einf462/ISR3D_UQ/ISR3D/kernel/absmc/build/equilibration3D $CONFIGNAME
/projects/0/einf462/ISR3D_UQ/ISR3D/kernel/absmc/build/deployment3D-stent $CONFIGNAME
