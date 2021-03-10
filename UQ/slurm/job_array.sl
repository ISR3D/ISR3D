#!/bin/bash

#### Job Array Setting
# Job Name and Files
#SBATCH -J SMCParaUQ

# Output and error
#SBATCH -N 1
#SBATCH --exclusive

# Wall time in format Day-hour:minutes:seconds
#SBATCH --time=0-00:30:00

## grant
#SBATCH --partition=normal

## Set job array,
#SBATCH --array=1-128%20



### Set number of samples for SA
# General parameter setting
NumSample=128
NumInput=4
SampleArray=("A")
Outputname=UQtest
### Actual UQ run setting

# Better to stay with N*num_sample
NumSimu=$SLURM_ARRAY_TASK_MAX
CurSimu=$(($SLURM_ARRAY_TASK_ID - 1))


### Pre-setting of Muscle3 and module loading
### The directory should now: Experiment_fold
###							: job_array.sh updatefunc.py printfunc.py
### 						: A 
###							: A_0
###							: input.ymmsl

### Individual Job pre-setting
ModValue=$((${CurSimu}%${NumSample}))

# Case start with A
MatValue=${SampleArray[0]}
MatValueSub=${SampleArray[0]}_${ModValue}

# Go to the corresponding inputs directory
cd ${Outputname}/${MatValue}/${MatValueSub}

#######################################################################
########################## Start Simulation ###########################
#######################################################################
echo $(date -u) "kill previous muscle"
killall muscle_manager


# Setting: activate C++&python muscle and load GCC/OpenMPI
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/muscle3/lib
source /home/dong/anaconda3/etc/profile.d/conda.sh
conda activate ISR3D


module load pre2019
module load OpenMPI/2.1.1-GCC-6.4.0-2.28

# Sleep for a while to make sure module loaded
sleep 1

# Run SMC and flow model
echo $(date -u) "Simulation started"
muscle_manager ./input_stage4.ymmsl &

sleep 1

# Start submodels​
export OMP_NUM_THREADS=24
/projects/0/einf462/ISR3D_UQ/ISR3D/build/smc --muscle-instance=smc                     >'smc.log' 2>&1  &
/projects/0/einf462/ISR3D_UQ/ISR3D/build/voxelizer --muscle-instance=voxelizer         >'vol.log' 2>&1  &
/projects/0/einf462/ISR3D_UQ/ISR3D/build/distributor --muscle-instance=distributor     >'dis.log' 2>&1  &
/projects/0/einf462/ISR3D_UQ/ISR3D/build/collector --muscle-instance=collector         >'col.log' 2>&1  &
srun -n 16 /projects/0/einf462/ISR3D_UQ/ISR3D/build/FlowController3D --muscle-instance=flow       >'flow.log' 2>&1

wait

echo $(date -u) "Whole simulation finished"

