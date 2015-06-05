#!/bin/bash



# for each directory Txxx
# generate the text energy file e.xvg from the nvt.edr



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

for dir in `ls --color=none -d T[0-9][0-9][0-9]`; do
  cd $dir
  echo "11 17 0" | $gmx energy -f nvt.edr -o e.xvg
  rm -f \#*
  cd ..
done
