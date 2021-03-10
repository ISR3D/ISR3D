gawk '{$2=$4; $4=-$5; $5=$2; $7=$4; print $0}' polimi_balloon.dat > polimi_balloon_rotated.dat
