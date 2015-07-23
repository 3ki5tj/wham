#!/bin/bash



# for each directory Txxxdbl
# generate the text energy file e.xvg from the nvt.edr



if [ -d $HOME/lwork/gmx ]; then
  export gmxdir=$HOME/lwork/gmx
else
  export gmxdir=$HOME/work/gmx
fi

export GMXLIB=$gmxdir/gromacs5.0/src/top

if [ -f $gmxdir/gromacs5.0/buildiccdbl/bin/gmx_d ]; then
  export gmx=$gmxdir/gromacs5.0/buildiccdbl/bin/gmx_d
elif [ -f $gmxdir/gromacs5.0/buildgccdbl/bin/gmx_d ]; then
  export gmx=$gmxdir/gromacs5.0/buildgccdbl/bin/gmx_d
else
  export gmx=$gmxdir/gromacs5.0/buildgccdbl32/bin/gmx_d
fi

for dir in `ls --color=none -d T[0-9][0-9][0-9]dbl`; do
  cd $dir
  echo "11 0" | $gmx energy -f nvt.edr -o e.xvg
  rm -f \#*
  cd ..
done
