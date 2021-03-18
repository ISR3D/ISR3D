import subprocess
lib_list = ['numpy','ymmsl','sobol_seq','csv','seaborn','zenodo_get']
for lib_name in lib_list:
	try:
		import lib_name
	except ImportError:
		if lib_name == 'csv':
			print(lib_name,' Module not installed')
			subprocess.run(['pip','install','python-csv']) 
		else:
			print(lib_name,' Module not installed')
			subprocess.run(['pip','install','%s'%lib_name])


import numpy as np
import ymmsl
import sobol_seq
import csv
import os
import seaborn as sns
import zenodo_get


# Transform the normalized sample matrix to ranges of uncertain parameters
def dim_transform(sobol_vector,uncertain_list):
    
    dim = len(uncertain_list)
    for num_dim in range(dim):
        para_max = uncertain_list[num_dim].get('max')
        para_min = uncertain_list[num_dim].get('min')
        sobol_vector[:,num_dim] = para_min + (para_max-para_min)*sobol_vector[:,num_dim]
        
    return sobol_vector

####################################################################################
##### Sample generation and UQ campaign creation (including instances folder)#######
####################################################################################
# Note:
#  This is used to generate UQ samples for only four biological parameters:
#	1) Endothelium endpoint 2)smc max stran 3)balloon extension 4) Fenestration probability

# Naming of Folder and files for samples
# Level 0:    UQtest                                    (UQ campaign name)
# Level 1:    UQtest/A                                  (sample matrix of sobol sequence)
# Level 2:    UQtest/A/A_X where X vary from 1 -> N     (N: number of samples)
# Level 3:    UQtest/A/A_X/input.ymmsl


### Main function
# Number of samples for UQ
# Note that ISR3D is a computationally intensive application.
# Running 128 instances would need some cluster resources
# You can start with a small number, 16 for instances.
NumSample = 128

# Template path to the ymmsl file (relative path from ISR3D/Result/UQtest/ to ISR3D/UQ/template/input_stage4.ymmsl)
input_path = '../../UQ/template/'
input_ymmsl_filename = 'input_stage4.ymmsl'

# Output directory for UQ campagin folder and name
output_path = './'
experiment_name = 'UQtest'


# Read in the data of template ymmsl file
with open(input_path+input_ymmsl_filename,'r') as f:
    ymmsl_data = ymmsl.load(f)

# Take out the unchanged model part and need-for-change settings part for ymmsl model
model = ymmsl_data.model
settings = ymmsl_data.settings


# Set uncertain parameters and its ranges as a list
ymmsl_uncertain_parameters = [
        {
            'name': 'smc.endo_endpoint',
            'min': 10.0,
            'max': 20.0
        },
        {
            'name': 'smc.balloon_extension',
            'min': 0.5,
            'max': 1.5
        },
        {
            'name': 'smc.smc_max_strain',
            'min': 1.2,
            'max': 1.8
        },
        {
            'name': 'smc.fenestration_probability',
            'min': 0.0,# Calculate the lumen volume from (lumen_area_of_each_slice*depth_of_slice)
            'max': 0.1
        }]

# Count the total uncertain input dimensions (here 4 parameters)
num_uncer_para = len(ymmsl_uncertain_parameters)
print('Number of uncertain parameter: '+str(num_uncer_para))

# Generate sobel sequence range (0,1), save the file and transform to (min,max)
A = sobol_seq.i4_sobol_generate(num_uncer_para,NumSample)
A = dim_transform(A,ymmsl_uncertain_parameters)
np.savetxt("A.csv",A)

# Create corresponding directory and folders
try:
    os.mkdir(output_path+experiment_name)
except OSError:
    print ("Creation of the directory %s failed" % output_path+experiment_name)
else:
    print ("Successfully created the directory %s" % output_path+experiment_name)


# A: Replace the corresponding value within the dict and output the file
os.mkdir(output_path+experiment_name+'/A')
checklist = ['A']

for n in range(NumSample):
    sample_path = output_path+experiment_name+'/A'+'/A_'+str(n)
    os.mkdir(sample_path)
    
    # Generate file for ymmsl 
    num_para = 0
    for para in ymmsl_uncertain_parameters:
        settings[para.get('name')] = float(A[n,num_para])
        num_para = num_para + 1
        
    config = ymmsl.Configuration(model, settings)
    with open(sample_path+'/input_stage4.ymmsl', 'w') as f:
        ymmsl.save(config, f)

print('ymmsl input for each UQ instance has been generated')

####################################################################################
##### Run shell script to broadcast other input files to each sample folder#########
####################################################################################
import subprocess
# Download Other input files from Zenodo
print('Start to download other input files for ISR3D from Zenodo')
subprocess.run(['wget https://zenodo.org/record/4603912/files/stage3.test_vessel.dat'],shell = True)
subprocess.run(['wget https://zenodo.org/record/4603912/files/stage3.test_vessel_nb.dat'],shell = True)
subprocess.run(['wget https://zenodo.org/record/4603912/files/test_vessel_centerline.csv'],shell = True)
print('Start to broadcast the input to each UQ instance directory')
# Template path to the ymmsl file (relative path from ISR3D/Result/UQtest/ to ISR3D/UQ/function/BCastStage3.sh)
pass_arg = str(NumSample)
subprocess.run(['bash','../../UQ/function/BCastStage3.sh', '%s'%pass_arg]) 
print('Sample generation done')