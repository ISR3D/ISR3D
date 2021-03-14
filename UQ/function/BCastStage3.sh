### Set number of samples for SA
# General parameter setting

echo $1
SampleArray=("A")
Outputname=UQtest
echo ${Outputname}
for Matrix in "${SampleArray[@]}"
do
	if [ $Matrix = A ]
	then
		for NSample in `seq 0 $(($1 - 1))`
		do
			cd ${Outputname}/${Matrix}/${Matrix}_${NSample}
			cp ../../../stage3.test_vessel.dat .
			cp ../../../stage3.test_vessel_nb.dat .
			cp ../../../test_vessel_centerline.csv .
			cd ../../../
		done
	fi
done
echo Broadcast finished



