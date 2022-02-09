import numpy as np
import matplotlib.pyplot as plt
import SimpleITK as sitk


a = np.load('../Dataset/numpyfile_of_binary_geometry.npy')
a[a!=0]=2
a[a==0]=1
a[a==2]=0
print(a.shape)

plt.imshow(a[:,:,45])
plt.colorbar()
plt.show()

# a = np.int_(a)

filtered_image = sitk.GetImageFromArray(a)

print(filtered_image)
sitk.WriteImage(filtered_image, "tmp/TestGeo.vtk")


# print(a.shape)
