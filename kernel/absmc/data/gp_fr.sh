#!/bin/bash

if [ $# -ne 1 ]
then
  echo "Usage: `basename $0` file"
  exit 1
fi  

FNAME=$1

gnuplot -persist << EOF

set terminal png enhanced size 1200,800
set output "$FNAME.png"

set title "physical solver: force residuals"

set mxtics
set mytics
set grid
set grid mxtics mytics

set logscale y
set xlabel "total time"
set ylabel "force residuals"

set x2label "total number of iterations"
set x2tics

plot "$FNAME" u 2:5 axes x1y1 t "maximum norm" w lp,\
"" u 2:13 axes x1y1 t "l2 norm" w lp;

EOF

exit 0


#"" u 1:8 axes x2y1 w l lw 2 t "eps*F0",\