#!/bin/bash

if [ $# -ne 2 ]
then
  echo "Usage: `basename $0` work_dir frame_rate"
  exit 1
fi  

OUT_DIR=$1
FRAME_RATE=$2
mencoder "mf://$OUT_DIR/*.png" -mf fps=$FRAME_RATE -o $OUT_DIR/output.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800  
