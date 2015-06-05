#!/bin/bash


# for each TxxxPxxx directory
# create an energy file ev.xvg from npt.edr



if [ -d $HOME/work/gmx ]; then
  export gmxdir=$HOME/work/gmx
else
  export gmxdir=$HOME/lwork/gmx
fi

export GMXLIB=$gmxdir/gromacs5.0/src/top

if [ -f $gmxdir/gromacs5.0/buildicc/bin/gmx ]; then
  export gmx=$gmxdir/gromacs5.0/buildicc/bin/gmx
elif [ -f $gmxdir/gromacs5.0/buildgcc/bin/gmx ]; then
  export gmx=$gmxdir/gromacs5.0/buildgcc/bin/gmx
else
  export gmx=$gmxdir/gromacs5.0/buildgcc32/bin/gmx
fi

for dir in `ls --color=none -d T[0-9][0-9][0-9]P*`; do
  cd $dir
  echo "11 16 21 0" | $gmx energy -f npt.edr -o ev.xvg
  rm -f \#*
  cd ..
done
