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

if [ "$1" -eq "2" ] ; then
  export mydirs="T310P0.95 T320P0.9 T330P0.85 T340P0.8 T350P0.75 T360P0.7 T370P0.65 T380P0.6 T390P0.55 T400P0.5 T410P0.45 T420P0.4"
elif [ "$1" -eq "3" ] ; then
  export mydirs="T310P2 T320P3 T330P4 T340P5 T350P6 T360P7 T370P8 T380P9 T390P10 T400P11 T410P12 T420P13"
else
  export mydirs=`ls --color=none -d T[0-9][0-9][0-9]P*`
fi

for dir in $mydirs; do
  cd $dir
  echo "11 16 21 0" | $gmx energy -f npt.edr -o ev.xvg
  rm -f \#*
  cd ..
done
