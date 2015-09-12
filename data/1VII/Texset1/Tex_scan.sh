#!/bin/bash



# for each directory Txxx
# generate the text energy file e.xvg from the nvt.edr



if [ -d $HOME/lwork/gmx/gromacs5.0 ]; then
  export gmxdir=$HOME/lwork/gmx
else
  export gmxdir=$HOME/work/gmx
fi

export GMXLIB=$gmxdir/gromacs5.0/src/top

if [ -f $gmxdir/gromacs5.0/buildicc/bin/gmx ]; then
  export gmx=$gmxdir/gromacs5.0/buildicc/bin/gmx
elif [ -f $gmxdir/gromacs5.0/buildgcc/bin/gmx ]; then
  export gmx=$gmxdir/gromacs5.0/buildgcc/bin/gmx
else
  export gmx=$gmxdir/gromacs5.0/buildgcc32/bin/gmx
fi

for dir in T450 T500 T550 T600 T650 T700 T750 T800 T850 T900 T950 T1000; do
  cd $dir
  echo "11 0" | $gmx energy -f nvt.edr -o e.xvg
  rm -f \#*
  cd ..
done
