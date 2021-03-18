#!/bin/bash

#### Job Array Setting
# Job Name and Files
#SBATCH -J SMCParaUQ

# Output and error
#SBATCH -N 1
#SBATCH --exclusive

# Wall time in format Day-hour:minutes:seconds
#SBATCH --time=0-00:30:00

## Partition, Important!!!: Need to be adapted your setting
#SBATCH --partition=normal

## Set job array,
#SBATCH --array=1-128%20


###################################################
########### Job array setting #####################
###################################################
# General parameter setting
NumSample=128
NumInput=4
SampleArray=("A")
Outputname=UQtest

# Fetch the total number of job, and current numbering of job
NumSimu=$SLURM_ARRAY_TASK_MAX
CurSimu=$(($SLURM_ARRAY_TASK_ID - 1))


### The directory should now: Experiment_fold
### 						: A 
###							: A_0
###							: input.ymmsl


# Find out the name of UQ instance
ModValue=$((${CurSimu}%${NumSample}))
MatValue=${SampleArray[0]}
MatValueSub=${SampleArray[0]}_${ModValue}

# Go to the corresponding inputs directory
cd ${Outputname}/${MatValue}/${MatValueSub}

#######################################################################
########################## Start Simulation ###########################
#######################################################################
echo $(date -u) "kill previous muscle if existing"
killall muscle_manager


# Setting: activate C++& python muscle  
# Important!!!: Need to be adapted your setting
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/muscle3/lib
source directory_to_anaconda3/anaconda3/etc/profile.d/conda.sh
conda activate Python_Env_Name

# Setting: load GCC/OpenMPI
# Important!!!: Need to be adapted your setting
module load OpenMPI_version

# Sleep for a while to make sure module loaded
sleep 1

# Run SMC and flow model
echo $(date -u) "Simulation started"
muscle_manager ./input_stage4.ymmsl &

sleep 1

# Start submodels​
# Important!!!: The directory of submodel needed to be adapted your setting!
# Important!!!: OMP_NUM_THREADS=XX adapted your setting!
# Important!!!: srun -n XX  adapted your setting!
export OMP_NUM_THREADS=XX
Directory_to_ISR3D/ISR3D/build/smc --muscle-instance=smc                     >'smc.log' 2>&1  &
Directory_to_ISR3D/ISR3D/build/voxelizer --muscle-instance=voxelizer         >'vol.log' 2>&1  &
Directory_to_ISR3D/ISR3D/build/distributor --muscle-instance=distributor     >'dis.log' 2>&1  &
Directory_to_ISR3D/ISR3D/build/collector --muscle-instance=collector         >'col.log' 2>&1  &
srun -n XX directory_to_ISR3D/ISR3D/build/FlowController3D --muscle-instance=flow       >'flow.log' 2>&1

wait

echo $(date -u) "Whole simulation finished"

