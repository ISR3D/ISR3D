### Set number of samples for SA
# General parameter setting
NumSample=128
SampleArray=("A")
Outputname=UQtest
echo ${Outputname}
for Matrix in "${SampleArray[@]}"
do
	if [ $Matrix = A ]
	then
		for NSample in `seq 0 $((${NumSample} - 1))`
		do
			cd ${Outputname}/${Matrix}/${Matrix}_${NSample}
			cp ../../../../../UQ/Input/stage3.test_vessel.dat .
			cp ../../../../../UQ/Input/stage3.test_vessel_nb.dat .
			cp ../../../../../UQ/Input/test_vessel_centerline.csv .
			cd ../../../
		done
	fi
done
echo Broadcast finished



