#!/bin/bash
#Usage: loadmuscle.sh configname.ymmsl
#NB: has to be run from the run directory two folders down

CONFIGNAME=$1
ISR3D_HOME=/home/pavel/src/ISR3D/
#cd ISR/ISR3D_MUSCLE3/run/LAD1785/
. ~/src/muscle3/muscle3_venv/bin/activate
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/pavel/apps/muscle3/lib

killall muscle_manager smc voxelizer distributor collector FlowController3D

#muscle_manager ../../cxa/isr_tinier_test_vessel.ymmsl
muscle_manager $ISR3D_HOME/cxa/$CONFIGNAME &

sleep 1

# Start submodels​

#openMP submodel
$ISR3D_HOME/build/smc --muscle-instance=smc                     >'smc_cout.log' 2>&1  &

#Helper modules to convert data between openMP and MPI parts
$ISR3D_HOME/build/voxelizer --muscle-instance=voxelizer         >'vol_cout.log' 2>&1  &
$ISR3D_HOME/build/distributor --muscle-instance=distributor     >'dis_cout.log' 2>&1  &
$ISR3D_HOME/build/collector --muscle-instance=collector         >'col_cout.log' 2>&1  &

#MPI submodel
mpirun -np 4 $ISR3D_HOME/build/FlowController3D --muscle-instance=flow       >'flow_cout.log' 2>&1

wait

killall muscle_manager smc voxelizer distributor collector FlowController3D
