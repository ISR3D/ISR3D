import csv
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# lattice dimensions
nx = 131
ny = 46
nz = 42

# lattice spacing
dx = 0.015

field1d = zeros (nx * ny * nz)

with open("./wssfield0000.dat", newline='') as csvfile:
    datareader = csv.reader(csvfile, delimiter=' ')
    for row in datareader:
        for index in range(nx * ny * nz):
            field1d[index] = float(row[index])

print("file read to 1d array")

field3d = field1d.reshape((nx,ny,nz))
print("1d array converted to 3d")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# x,y,z = meshgrid(range(nx), range(ny), range(nz))
# ax.scatter(x,y,z, c=field3d.flat)

x,y,z = meshgrid((40), range(ny), range(nz))
ax.scatter(x,y,z, c=field3d[40].flat)
plt.show()
