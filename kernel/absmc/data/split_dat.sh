#!/bin/bash

# shell script to split absmc dat files containing agents of several types
# into multiple files containing only one type of agent each.
# multiple input files allowed.

if [ $# -lt 1 ]
then
  echo "Usage: `basename $0` FILES. Multiple input files allowed."
  exit 1
fi  

for FN in $@
do
    echo "Processing $FN ..."
    
    # output file names
    BASENAME=${FN%".dat"}
    FN_SMC3D=$BASENAME"_SMC3D.dat"
    FN_IEL3D=$BASENAME"_IEL3D.dat"
    
    # typeId constants as defined in src/core/agentTypeId.h
    tAny=0
    tSMC2D=1
    tIEL2D=2
    tObstacle2D=3
    tSMC3D=4
    tIEL3D=5
    tObstacle3D=6
    
    awk '$1=="'$tSMC3D'" {print}' $FN > $FN_SMC3D
    awk '$1=="'$tIEL3D'" {print}' $FN > $FN_IEL3D
done

exit 0
