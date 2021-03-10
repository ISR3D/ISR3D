#these commands make a .dat from .csv for stents and other obstacles
#use csvArteryLoader3D for arteries

gawk -F "," '{print 6 " " $1 " " $2 " " $3 " " $1 " " $2 " " $3 " 0.015 0 0 0"}' RCA929-Job-v7r5_CU28_stent.csv > RCA929_v2.dat
gawk -F "\t" '{print 6 " " $1 " " $2 " " $3 " " $1 " " $2 " " $3 " 0.03 0 0 0"}' balloonagents.csv > polimi_balloon.dat
gawk -F "\t" '{print 6 " " $1 " " $2 " " $3 " " $1 " " $2 " " $3 " 0.03 0 0 0"}' synergy-LAD-1785.csv > synergy-LAD-1785.dat
