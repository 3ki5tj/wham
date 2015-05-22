#!/bin/bash


# for each TxxxPxxx directory
# create an energy file ev.xvg from npt.edr



export gmxdir=$HOME/work/gmx
export GMXLIB=$gmxdir/gromacs5.0/src/top

for dir in `ls --color=none -d T[0-9][0-9][0-9]P*`; do
  cd $dir
  echo "11 16 21 0" | $gmxdir/gromacs5.0/buildicc/bin/gmx energy -f npt.edr -o ev.xvg
  rm -f \#*
  cd ..
done
