import numpy as np
import ymmsl
import yaml
import sobol_seq
import os

# Note:
#  This is used to generate UQ samples for only four biological parameters:
#	1) Endothelium endpoint 2)smc max stran 3)Blood flow velocity 4) Fenestration probability
#  All the parameters in cfg file are disabled (Take average value)
#  Please use the muscle3 python environment 0.3.2 for execution

# Quasi-Monte Carlo method (Sobol seqence) is applied to generate samples here.
# The structures of directory and sample names are implemented in a Saltelli's way which would be easier for
# Sobol sensitivity analysis (without surrogate model, pure MC implementation). If you just want generate
# training samples for UQ with qMC, disable the  B/AB/BA folder generation.

# Naming of Folder and files for samples
# Level 1:    A, B, AB, B
# Level 2:    A/B: 1 -> N       AB/BA: 1 -> D   (N:SA sample, D: inputs dimensions)
# Level 3:    A/B: cfg,ymmsl    AB/BA: 1 -> N
# Level 4:                      AB/BA: cfg,ymmsl

# e.g: dircatory of inputs files for a sample in sample matrix A:
#      A/A_1/input_stage1.cfg, A/A_1/input_stage4.ymmsl
# e.g: dircatory of inputs files for a sample in sample matrix AB:
#      AB/AB_1/AB_1_3/input_stage1.cfg, AB/AB_1/AB_1_3/input_stage4.ymmsl





# Read_ymmsl_file, TODO may remove the yaml reading/writing
def read_ymmsl_file(path,filename):
	final_path = path + filename
	with open(final_path) as f:
		ymmsl_data = ymmsl.load(f)
		
	with open(final_path) as f:
		yaml_data = yaml.load(f,Loader=yaml.FullLoader)

	return ymmsl_data,yaml_data

# Read cfg file
def read_cfg_file(path,filename):
	final_path = path + filename
	with open(final_path) as f:
		cfg_data = f.readlines() 

	return cfg_data

# Transfrom the cfg plain text to list and then dict
def clean_cfg_data(data):
	
	datalines=[]
	for line in data:
		#remove linebreaks
		newline = line.replace('\n','')

		#find start of comments
		location = newline.find('#')
		if location >= 0:
			newline = newline[0:location].strip()

		if newline is not '':
			datalines.append(newline)
	
	result_dict={}
	for line in datalines:
		location = line.find(' ')
		if location >= 0:
			name_vara = line[0:location].strip()
			value_vara = line[location:].strip()

		result_dict[name_vara] = value_vara

	return result_dict

# Transfrom dict data to cfg data (plain text with space)
def from_dict_to_cfg(some_dict):
	str_list=[]
	for key in some_dict:
		#print(key)
		#print(some_dict[key])
		str_list.append(key+' '+str(some_dict[key]))

	return str_list

def write_cfg_file(path,filename,config_data):
	
	f=open(path+filename, 'w')
	for str_data in config_data:
		f.write(str_data+'\n')
	f.close()

# save ymmsl file
def write_ymmsl_file(path,filename,config_data):
	with open(path+filename, 'w') as f:
		ymmsl.save(config_data, f)

# Generate Sobel sequence for both A and B
def sobol_gen(dimensions,num_sample):
	gen_vector = sobol_seq.i4_sobol_generate(dimensions*2,num_sample)
	return gen_vector[:,:dimensions],gen_vector[:,dimensions:]

# Transfrom the sobel sequence with min/max in dict of uncertain parameters
def dim_transform(sobel_vector,uncertain_list_1,uncertain_list_2):

	dim_1 = len(uncertain_list_1)
	dim_2 = len(uncertain_list_2)
	dim = dim_1 + dim_2

	for num_dim in range(dim):
		
		if num_dim < dim_1:
			para_max = uncertain_list_1[num_dim].get('max')
			para_min = uncertain_list_1[num_dim].get('min')
			sobel_vector[:,num_dim] = para_min + (para_max-para_min)*sobel_vector[:,num_dim]
		else:
			para_max = uncertain_list_2[num_dim-dim_1].get('max')
			para_min = uncertain_list_2[num_dim-dim_1].get('min')
			sobel_vector[:,num_dim] = para_min + (para_max-para_min)*sobel_vector[:,num_dim]

	return sobel_vector

# This function manually fix the:
#      1)  outer_radius = tunica_width + inner_radius
#      2)  deploy_depth = deploy_depth * inner radius
#      3)  curvature_radius = 1/curvature * LX
# take whole matrix A,B or A_B, B_A
# This is seperated from XX_fix function because its in between the uncertain input list
def special_para_tran(sobel_vector):

	sobel_vector[:,2] = sobel_vector[:,2] + sobel_vector[:,1]
	sobel_vector[:,4] = sobel_vector[:,4] * sobel_vector[:,1]
	sobel_vector[:,3] = (1/sobel_vector[:,3]) * sobel_vector[:,0] * stent2Lx

	return sobel_vector

# Replace the name of tunica_width to outer)radius for output, Be cautious!!!!!
def name_replace(list_dict):

	list_dict[2]['name'] = 'outer_radius'

	return list_dict


# Fix the problem of LX = 2*stent_length
def lx_fix(string_key,target_dict,target_matrix,num):
	target_dict[string_key] = target_matrix[num,0]*stent2Lx
	return target_dict

# Fix the problem of helix_ptich = inner_radius * helix_const	
def helix_fix(string_key,target_dict,target_matrix,num):
	target_dict[string_key] = target_matrix[num,1]*helix_const
	return target_dict

# Fix the problem of num_strut, integer from 6 to 18 as its linear to the inner radius from 0.5 to 1.5
def strut_fix(string_key,target_dict,target_matrix,num):
	radius_temp = target_matrix[num,1]
	strut_app = ((radius_temp - inner_radius_min)/(inner_radius_max - inner_radius_min))*(strut_max-strut_min)+strut_min

	target_dict[string_key] = int(strut_app)
	return target_dict



#########  Main function ##########
# Basic settings
N_sample = 512
input_path = './'

input_ymmsl_filename = 'input_stage4.ymmsl'
input_cfg_filename = 'input_stage1.cfg'

output_path = './'
experiment_name = 'UQ'

# Read existing sample data (disabled)
read_sample_data = 0

#matrix activation (for data generation, only A is needed.)
acti_A = 1
acti_B = 0
acti_AB = 0
acti_BA = 0

###### Const for setting
# Data comes from test file: tinier_test_vessl const = 2.0/0.4
helix_const = 12

inner_radius_max = 1.5
inner_radius_min = 0.5

strut_max = 18
strut_min = 6

stent2Lx = 1.5

# Read Data from ymmsl and cfg file, and process cfg data to dict
print('Start to read data from cfg/ymmsl template')
ymmsl_data,yaml_data = read_ymmsl_file(input_path,input_ymmsl_filename)

cfg_data = read_cfg_file(input_path,input_cfg_filename)
cfg_data = clean_cfg_data(cfg_data)

print('Data reading finished')
# Take out the unchanged model part and need-for-change settings part for ymmsl model
model = ymmsl_data.model
settings = ymmsl_data.settings


# Manually set the parameter 
# (the sequence will affect its position in the vector of sample)
print('Uncertain parameter reading')

# All the CFG parameter here will be set to average
cfg_uncertain_parameters = [
		{
			# To be confirmed its asssociation with lx,spiral_step,num_strut (checked)
			'name': 'stent_length',
			'min': 8.0,
			'max': 16.0
		},
		{
			# inner_radius = lumen_width/2 (checked)
			'name': 'inner_radius',
			'min': inner_radius_min,
			'max': inner_radius_max
		},
		{
			# outer_radius = inner_radius + tunica_width (checked)
			'name': 'tunica_width',
			'min': 0.2,
			'max': 0.5
		},
		{
			# Note that the parameter ranges here is actually curvature
			# The name curvature_radius is used just for easy generation of input files
			'name': 'curvature_radius',
			'min': 0.0,
			'max': np.pi/2
		},
		{
			# To be confirmed, deployment_depth(in percentage)*inner_radius (checked)
			'name': 'deployment_depth',
			'min': 0.05,
			'max': 0.2
		}]


ymmsl_uncertain_parameters = [
        {
            'name': 'flow.flow_velocity',
            'min': 0.133,
            'max': 0.399
        },
        {
        # (checked)
            'name': 'smc.endo_endpoint',
            'min': 10.0,
            'max': 20.0
        },
        # {
        # 	# To be confirmed the value
        #     'name': 'smc.balloon_extension',
        #     'min': 0.5,
        #     'max': 1.5
        # }, 
        # {
        #     'name': 'smc.smc_mean_rad',
        #     'min': 0.012,
        #     'max': 0.018
        # },
        {
            'name': 'smc.smc_max_strain',
            'min': 0.446,
            'max': 0.785
        },
        {
            'name': 'smc.fenestration_probability',
            'min': 0.02,
            'max': 0.10
        }]



# Count the total uncertain input dimensions (here only 4 parameter, cfg file is excluded)
num_uncer_para = len(ymmsl_uncertain_parameters)# + len(cfg_uncertain_parameters)
print('Number of uncertain parameter: '+str(num_uncer_para))


# Generate sobel sequence range (0,1) and transform to (min,max)
if read_sample_data == 0:
	A,B = sobol_gen(num_uncer_para,N_sample)
	np.savetxt("A.csv",A)
	np.savetxt("B.csv",B)
else:
	A = loadtxt("A.csv")
	B = loadtxt("B.csv")

Acomp = np.ones((N_sample,len(cfg_uncertain_parameters)))/2 # Generate compensation matrix
Bcomp = np.ones((N_sample,len(cfg_uncertain_parameters)))/2

A = np.concatenate((Acomp,A),axis=1)
B = np.concatenate((Bcomp,B),axis=1)

A = dim_transform(A,cfg_uncertain_parameters,ymmsl_uncertain_parameters)
B = dim_transform(B,cfg_uncertain_parameters,ymmsl_uncertain_parameters)

# Clean Special Relation in uncertain inputs and practical inputs
A = special_para_tran(A)
B = special_para_tran(B)

cfg_uncertain_parameters = name_replace(cfg_uncertain_parameters)
print('Sobel Sequence generation and data processing done')


# Create corresponding directory and folders
try:
    os.mkdir(output_path+experiment_name)
except OSError:
    print ("Creation of the directory %s failed" % output_path+experiment_name)
else:
    print ("Successfully created the directory %s" % output_path+experiment_name)



# A: Replace the corresponding value within the dict and output the file

if acti_A == 1:
	os.mkdir(output_path+experiment_name+'/A')
	checklist = ['A']

	for n in range(N_sample):
		sample_path = output_path+experiment_name+'/A'+'/A_'+str(n)
		os.mkdir(sample_path)

		# Generate file for cfg
		num_para = 0
		for para in cfg_uncertain_parameters:
			cfg_data[para.get('name')] = A[n,num_para]
			num_para = num_para + 1

		cfg_data = lx_fix('lX',cfg_data,A,n)
		cfg_data = helix_fix('helix_pitch',cfg_data,A,n)
		cfg_data = strut_fix('num_struts',cfg_data,A,n)
		cfg_data['stent_curvature_radius'] = cfg_data['curvature_radius']

		# This fix where boundary_trim = 2*outer_radius*sin(Lx/2/curvature_radius)
		local_boundary_trim = 2*cfg_data['outer_radius']*np.sin(cfg_data['lX']/2/cfg_data['curvature_radius'])

		output_str = from_dict_to_cfg(cfg_data)
		write_cfg_file(sample_path+'/', 'input_stage1.cfg',output_str)

		# Generate file for ymmsl 
		for para in ymmsl_uncertain_parameters:
			settings[para.get('name')] = float(A[n,num_para])
			num_para = num_para + 1

		settings['smc.boundary_trim_distance'] = float(local_boundary_trim)

		config = ymmsl.Configuration(model, settings)
		write_ymmsl_file(sample_path+'/','input_stage4.ymmsl',config)

	print('Wrting and saving procedure done: A')

# B: Replace the corresponding value within the dict and output the file
if acti_B == 1:
	os.mkdir(output_path+experiment_name+'/B')
	checklist.append('B')

	for n in range(N_sample):
		sample_path = output_path+experiment_name+'/B'+'/B_'+str(n)
		os.mkdir(sample_path)

		# Generate file for cfg
		num_para = 0
		for para in cfg_uncertain_parameters:
			cfg_data[para.get('name')] = B[n,num_para]
			num_para = num_para + 1

		cfg_data = lx_fix('lX',cfg_data,B,n)
		cfg_data = helix_fix('helix_pitch',cfg_data,B,n)
		cfg_data = strut_fix('num_struts',cfg_data,B,n)
		cfg_data['stent_curvature_radius'] = cfg_data['curvature_radius']
		local_boundary_trim = 2*cfg_data['outer_radius']*np.sin(cfg_data['lX']/2/cfg_data['curvature_radius'])

		output_str = from_dict_to_cfg(cfg_data)
		write_cfg_file(sample_path+'/', 'input_stage1.cfg',output_str)

		# Generate file for ymmsl 
		for para in ymmsl_uncertain_parameters:
			settings[para.get('name')] = float(B[n,num_para])
			num_para = num_para + 1

		settings['smc.boundary_trim_distance'] = float(local_boundary_trim)

		config = ymmsl.Configuration(model, settings)
		write_ymmsl_file(sample_path+'/','input_stage4.ymmsl',config)

	print('Wrting and saving procedure done: B')

# AB: Replace the corresponding value within the dict and output the file
if acti_AB == 1:
	os.mkdir(output_path+experiment_name+'/AB')
	checklist.append('AB')

	for d in range(num_uncer_para):
		dim_path = output_path+experiment_name+'/AB'+'/AB_'+str(d)
		os.mkdir(dim_path)

		AB = np.copy(A)
		AB[:,d+len(cfg_uncertain_parameters)] = np.copy(B[:,d+len(cfg_uncertain_parameters)])

		for n in range(N_sample):
			sample_path = dim_path + '/AB_'+str(d) + '_' + str(n)
			os.mkdir(sample_path)

			# Generate file for cfg
			num_para = 0
			for para in cfg_uncertain_parameters:
				cfg_data[para.get('name')] = AB[n,num_para]
				num_para = num_para + 1

			cfg_data = lx_fix('lX',cfg_data,AB,n)
			cfg_data = helix_fix('helix_pitch',cfg_data,AB,n)
			cfg_data = strut_fix('num_struts',cfg_data,AB,n)
			cfg_data['stent_curvature_radius'] = cfg_data['curvature_radius']
			local_boundary_trim = 2*cfg_data['outer_radius']*np.sin(cfg_data['lX']/2/cfg_data['curvature_radius'])

			output_str = from_dict_to_cfg(cfg_data)
			write_cfg_file(sample_path+'/', 'input_stage1.cfg',output_str)

			# Generate file for ymmsl 
			for para in ymmsl_uncertain_parameters:
				settings[para.get('name')] = float(AB[n,num_para])
				num_para = num_para + 1

			settings['smc.boundary_trim_distance'] = float(local_boundary_trim)

			config = ymmsl.Configuration(model, settings)
			write_ymmsl_file(sample_path+'/', 'input_stage4.ymmsl',config)

	print('Wrting and saving procedure done: AB')

# BA: Replace the corresponding value within the dict and output the file
if acti_BA == 1:
	os.mkdir(output_path+experiment_name+'/BA')
	checklist.append('BA')

	for d in range(num_uncer_para):
		dim_path = output_path+experiment_name+'/BA'+'/BA_'+str(d)
		os.mkdir(dim_path)

		BA = np.copy(B)
		BA[:,d+len(cfg_uncertain_parameters)] = np.copy(A[:,d+len(cfg_uncertain_parameters)])

		for n in range(N_sample):
			sample_path = dim_path + '/BA_'+str(d) + '_' + str(n)
			os.mkdir(sample_path)
			
			# Generate file for cfg
			num_para = 0
			for para in cfg_uncertain_parameters:
				cfg_data[para.get('name')] = BA[n,num_para]
				num_para = num_para + 1

			cfg_data = lx_fix('lX',cfg_data,BA,n)
			cfg_data = helix_fix('helix_pitch',cfg_data,BA,n)
			cfg_data = strut_fix('num_struts',cfg_data,BA,n)
			cfg_data['stent_curvature_radius'] = cfg_data['curvature_radius']
			local_boundary_trim = 2*cfg_data['outer_radius']*np.sin(cfg_data['lX']/2/cfg_data['curvature_radius'])

			output_str = from_dict_to_cfg(cfg_data)
			write_cfg_file(sample_path+'/', 'input_stage1.cfg',output_str)

			# Generate file for ymmsl 
			for para in ymmsl_uncertain_parameters:
				settings[para.get('name')] = float(BA[n,num_para])
				num_para = num_para + 1

			settings['smc.boundary_trim_distance'] = float(local_boundary_trim)

			config = ymmsl.Configuration(model, settings)
			write_ymmsl_file(sample_path+'/', 'input_stage4.ymmsl',config)

	print('Wrting and saving procedure done: BA')


# Write out a ymmsl file as a checklsit to record the simulation
P={}
P_plus={}
G={}


for value in checklist:
	if value == 'A' or value == 'B':
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



