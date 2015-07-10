#!/bin/bash

export remote=oo1@lonestar.tacc.utexas.edu:scratch/wham/data/1VII/
#python prep.py -T 300 -P 1.0 --nt=4

mydirs="T320P3 T330P4 T340P5 T350P6 T360P7 T370P8 T380P9 T390P10 T400P11 T410P12 T420P13 T430P14"

for dir in $mydirs; do
  #echo $dir
  #ls $dir/*.xvg
  rsync -avz $dir $remote
done


