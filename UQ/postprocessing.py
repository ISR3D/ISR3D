import subprocess
lib_list = ['numpy','csv','seaborn','matplotlib']
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
import csv
import os 
import seaborn as sns
import matplotlib.pyplot as plt

################################################################################
###### Fetch the result data and plot out the PDF and mean+/-STD figures #######
################################################################################
# Function of reading data
def LumenRead(path,numfile):
    
    resultlist = np.empty((0,39), float)
    for i in range(numfile):
        filename = os.path.join(path,('lumen_area_000'+str("{0:0=3d}".format(i))+'.csv'))
        reader = csv.reader(open(filename, "r"), delimiter='\t')
        x = list(reader)
        result = np.array(x[0][:-1]).astype("float")
        resultlist = np.append(resultlist, np.expand_dims(result,axis=0), axis=0)
        
    # print(resultlist.shape)
    return resultlist





# Set Directory and number of instance in the UQ campaign
data_root = './UQtest/A/'
time_step = 361

# Read the all subdirectory of UQ instances in a list
data_list = [os.path.join(data_root, item) for item in sorted(os.listdir(data_root))]

#Create an empty list and fetch the data in an loop
lumen_list = np.empty((0,time_step,39), float)
for item in data_list:
    print('Processing:',item) 
    Data = LumenRead(item,time_step)
    lumen_list = np.append(lumen_list,np.expand_dims(Data,axis=0),axis=0)

print(lumen_list.shape)	
np.save('LumenData',lumen_list)


# Calculate the lumen volume from (lumen_area_of_each_slice*depth_of_slice)
lumen_list = np.load('LumenData.npy')
lumen_vol = np.sum(lumen_list[:,:,:],axis=2) * 0.03125
fig = plt.figure()

# plt.plot(np.ones(128)*3.56055,np.linspace(0.0,2.5,128),label='0 days',c='k')
sns.kdeplot(lumen_vol[:,72],label='3 days',fill=True)
sns.kdeplot(lumen_vol[:,144],label='6 days',fill=True)
sns.kdeplot(lumen_vol[:,216],label='9 days',fill=True)
sns.kdeplot(lumen_vol[:,288],label='12 days',fill=True)
sns.kdeplot(lumen_vol[:,360],label='15 days',fill=True)
plt.legend(fontsize=10,loc=2)
plt.xlabel('Lumen volume of blood vessel ($mm^3$)',fontsize=12)
plt.ylabel('Probability density function',fontsize=12)
plt.savefig('./'+'pdf.png')
plt.show()



# plot mean+/-STD
days = np.zeros(361)
for i in range(361):
    days[i] = i/24

plt.plot(days, np.mean(lumen_vol,axis=0),label='Mean')
plt.plot(days, np.mean(lumen_vol,axis=0) - np.std(lumen_vol,axis=0), '--',color='r',label='Mean+/-STD')
plt.plot(days, np.mean(lumen_vol,axis=0) + np.std(lumen_vol,axis=0), '--',color='r')
plt.xlabel('Time (days)')
plt.ylabel('Lumen volume of blood vessel ($mm^3$)')
plt.legend()
plt.savefig('./'+'mean.png')
plt.show()