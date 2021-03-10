import numpy as np
import ymmsl
import yaml
import sobol_seq
import os

# Note:
#  This is used to generate UQ samples for only four biological parameters:
#	1) Endothelium endpoint 2)smc max stran 3)balloon extension 4) Fenestration probability

# Naming of Folder and files for samples
# Level 1:    A     (sample matrix of sobol sequence)
# Level 2:    A/A_X where X vary from 1 -> N     (N: number of samples)
# Level 3:    A/A_X/input.ymmsl

# e.g: dircatory of inputs files for a sample in sample matrix A: 
#     A/A_1/input_stage4.ysmmsl



# Read_ymmsl_file, TODO may remove the yaml reading/writing
def read_ymmsl_file(path,filename):
	final_path = path + filename
	with open(final_path) as f:
		ymmsl_data = ymmsl.load(f)
		
	with open(final_path) as f:
		yaml_data = yaml.load(f,Loader=yaml.FullLoader)

	return ymmsl_data,yaml_data


# save ymmsl file
def write_ymmsl_file(path,filename,config_data):
	with open(path+filename, 'w') as f:
		ymmsl.save(config_data, f)

# Generate Sobel sequence for both A and B
def sobol_gen(dimensions,num_sample):
	gen_vector = sobol_seq.i4_sobol_generate(dimensions,num_sample)
	return gen_vector

# Transfrom the sobel sequence with min/max in dict of uncertain parameters
def dim_transform(sobel_vector,uncertain_list):

	dim = len(uncertain_list)
	for num_dim in range(dim):
		para_max = uncertain_list[num_dim].get('max')
		para_min = uncertain_list[num_dim].get('min')
		sobel_vector[:,num_dim] = para_min + (para_max-para_min)*sobel_vector[:,num_dim]

	return sobel_vector



#########  Main function ##########
# Basic settings
N_sample = 128
input_path = './'

input_ymmsl_filename = 'input_stage4.ymmsl'

output_path = './'
experiment_name = 'UQtest'


# Read Data from ymmsl and cfg file, and process cfg data to dict
print('Start to read data from cfg/ymmsl template')
ymmsl_data,yaml_data = read_ymmsl_file(input_path,input_ymmsl_filename)
print('Data reading finished')

# Take out the unchanged model part and need-for-change settings part for ymmsl model
model = ymmsl_data.model
settings = ymmsl_data.settings

# Manually set the parameter 
# (the sequence will affect its position in the vector of sample)
print('Uncertain parameter reading')
ymmsl_uncertain_parameters = [
        {
        	# (checked)
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
            'min': 0.0,
            'max': 0.1
        }]



# Count the total uncertain input dimensions (here 4 parameters)
num_uncer_para = len(ymmsl_uncertain_parameters)
print('Number of uncertain parameter: '+str(num_uncer_para))


# Generate sobel sequence range (0,1), save the file and transform to (min,max)
A= sobol_gen(num_uncer_para,N_sample)
np.savetxt("A.csv",A)
A = dim_transform(A,ymmsl_uncertain_parameters)
print('Sobel Sequence generation and data processing done')


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

for n in range(N_sample):
	sample_path = output_path+experiment_name+'/A'+'/A_'+str(n)
	os.mkdir(sample_path)

	# Generate file for ymmsl 
	num_para = 0
	for para in ymmsl_uncertain_parameters:
		settings[para.get('name')] = float(A[n,num_para])
		num_para = num_para + 1

	config = ymmsl.Configuration(model, settings)
	write_ymmsl_file(sample_path+'/','input_stage4.ymmsl',config)

print('Wrting and saving procedure done: A')


# Write out a ymmsl file as a checklsit to record the simulation
P={}
P_plus={}
G={}


for value in checklist:
	if value == 'A':
		for i in range(N_sample):
			P[value+'_'+str(i)] = 0
		G[value] = P
		P={}
	else:
		for j in range(num_uncer_para):
			for i in range(N_sample):
				P[value+'_'+str(j)+'_'+str(i)] = 0
			P_plus[value+'_'+str(j)] = P
			P={}
		G[value] = P_plus
		P_plus={}

with open(output_path+experiment_name+'/'+'checklist.yaml', 'w') as f:
    data = yaml.dump(G, f)


 

