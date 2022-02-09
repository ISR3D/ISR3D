NumSample=128

for NSample in `seq 0 $((${NumSample} - 1))`
do
	echo ${NSample}
	### Preprocessing
	python binary2vtk.py ${NSample}
	#sleep 2s

	### Apply levelset and segmentation
	vmtkimageinitialization -ifile day10/binary/${NSample}.vtk -interactive 0 -method threshold -upperthreshold 1.1 -lowerthreshold 0.9 --pipe vmtklevelsetsegmentation -ifile day10/binary/${NSample}.vtk -iterations 300 -ofile day10/levelset/${NSample}.vti
	sleep 2s
	vmtkmarchingcubes -ifile day10/levelset/${NSample}.vti -ofile day10/model/${NSample}.vtp
	sleep 2s
	vmtksurfacesmoothing -ifile day10/model/${NSample}.vtp -passband 0.1 -iterations 30 -ofile day10/model_sm/${NSample}.vtp
	### Computer cross-sectional area
	vmtkcenterlineresampling -ifile model_centerline_sm.vtp -length 1.0 --pipe vmtkcenterlinesections -ifile day10/model_sm/${NSample}.vtp -ofile day10/area/${NSample}.vtp
	
	### Convert with paraview to csv file
	timeout 10 ../../../../Software/ParaView-5.8.1-MPI-Linux-Python3.7-64bit/bin/pvpython vtk2numpy.py ${NSample}
	
done

