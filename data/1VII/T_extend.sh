#!/bin/bash



# extend simulation time
# this script only needs to run for once



export gmxdir=$HOME/work/gmx
export GMXLIB=$gmxdir/gromacs5.0/src/top

for dir in `ls --color=none -d T[0-9][0-9][0-9]`; do
  cd $dir
  cp nvt.tpr nvt_old.tpr
  $gmxdir/gromacs5.0/buildicc/bin/gmx convert-tpr -s nvt.tpr -o nvt_new.tpr -extend 1980000
  mv nvt_new.tpr nvt.tpr
  rm -f \#*
  cd ..
done
