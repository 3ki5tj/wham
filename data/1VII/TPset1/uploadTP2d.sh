#!/bin/bash

export remote=oo1@lonestar.tacc.utexas.edu:scratch/wham/data/1VII/
#python prep.py -T 300 -P 1.0 --nt=4

mydirs="T300P1.1 T300P1.2 T300P1.3 T300P1.4 T300P1.5 T300P1.6 T300P1.7 T300P1.8 T300P1.9 T300P2 T300P2.1"
mydirs="$mydirs T310P1 T310P1.2 T310P1.3 T310P1.4 T310P1.5 T310P1.6 T310P1.7 T310P1.8 T310P1.9 T310P2 T310P2.1"
mydirs="$mydirs T320P1 T320P1.1 T320P1.3 T320P1.4 T320P1.5 T320P1.6 T320P1.7 T320P1.8 T320P1.9 T320P2 T320P2.1"
mydirs="$mydirs T330P1 T330P1.1 T330P1.2 T330P1.4 T330P1.5 T330P1.6 T330P1.7 T330P1.8 T330P1.9 T330P2 T330P2.1"
mydirs="$mydirs T340P1 T340P1.1 T340P1.2 T340P1.3 T340P1.5 T340P1.6 T340P1.7 T340P1.8 T340P1.9 T340P2 T340P2.1"
mydirs="$mydirs T350P1 T350P1.1 T350P1.2 T350P1.3 T350P1.4 T350P1.6 T350P1.7 T350P1.8 T350P1.9 T350P2 T350P2.1"
mydirs="$mydirs T360P1 T360P1.1 T360P1.2 T360P1.3 T360P1.4 T360P1.5 T360P1.7 T360P1.8 T360P1.9 T360P2 T360P2.1"
mydirs="$mydirs T370P1 T370P1.1 T370P1.2 T370P1.3 T370P1.4 T370P1.5 T370P1.6 T370P1.8 T370P1.9 T370P2 T370P2.1"
mydirs="$mydirs T380P1 T380P1.1 T380P1.2 T380P1.3 T380P1.4 T380P1.5 T380P1.6 T380P1.7 T380P1.9 T380P2 T380P2.1"
mydirs="$mydirs T390P1 T390P1.1 T390P1.2 T390P1.3 T390P1.4 T390P1.5 T390P1.6 T390P1.7 T390P1.8 T390P2 T390P2.1"
mydirs="$mydirs T400P1 T400P1.1 T400P1.2 T400P1.3 T400P1.4 T400P1.5 T400P1.6 T400P1.7 T400P1.8 T400P1.9 T400P2.1"
mydirs="$mydirs T410P1 T410P1.1 T410P1.2 T410P1.3 T410P1.4 T410P1.5 T410P1.6 T410P1.7 T410P1.8 T410P1.9 T410P2"

for dir in $mydirs; do
  #echo $dir
  #ls $dir/*.xvg
  rsync -avz $dir $remote
done


