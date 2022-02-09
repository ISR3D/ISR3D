import numpy as np
import os
import csv
import ymmsl
import pandas as pd

nx = 526
ny = 148
nz = 100

def load_solodata(filepath):

	wssfield = np.zeros(nx * ny * nz)

	with open(filepath, newline='') as csvfile:
		datareader = csv.reader(csvfile, delimiter=' ')
		for row in datareader:
			for index in range(nx * ny * nz):
				# print(index)
				wssfield[index] = float(row[index])

	wssfield = wssfield.reshape((nx,ny,nz))
	return wssfield



def LumenRead(path):

	filename = path+'wssfield0720.dat'
	reader = load_solodata(filename)

	return reader

def InputRead(path):

	filename = os.path.join(path,('input_stage4.ymmsl'))
	with open(filename) as f:
		ymmsl_data = ymmsl.load(f)

	settings = ymmsl_data.settings
	return settings





data_root = 'directory_to_wss_dat_file_generated_from_ISR3D'
time_step = 721

lumen_list = np.empty((0,526,142,90), float)
DataPath = data_root + '/'
Data = LumenRead(DataPath)
#lumen_list = np.append(lumen_list,np.expand_dims(Data,axis=0),axis=0)

print(Data.shape)	
np.save('NM_LastLumen_Max',Data)



