gawk '{$3=$4; $4=-$6; $6=$3; $7=$4; print $0}' RCA929-v7r5.dat > RCA929-rotated.dat