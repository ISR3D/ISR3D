#!/bin/bash

if [ $# -ne 1 ]
then
  echo "Usage: `basename $0` dir"
  exit 1
fi  

DIRNAME=$1

gnuplot -persist << EOF

set terminal png enhanced size 1200,800
set output "$DIRNAME/eq.png"

set mxtics
set mytics
set grid
set grid mxtics mytics

set logscale y
set xlabel "total time"
set ylabel "force residuals"

set x2label "total number of iterations"
set x2tics

plot "$DIRNAME/fr.dat" u 2:5 axes x1y1 t "max(force res)" w lp, "" u 2:13 axes x1y1 t "l2(force res)" w lp,\
"$DIRNAME/l2diff.dat" u 3:4 axes x1y1 t "l2(diff)" w lp, "" u 3:5:6 axes x1y1 t "nb dist" w yerr;

EOF

exit 0
