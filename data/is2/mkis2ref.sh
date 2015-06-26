#!/bin/sh

# make the data files for is2ref.gp

pdir=../../prog/is2
prog=$pdir/is2wham
histref=histref.dat

make -C $pdir $prog

if [ -f $histref ]; then
  opt=-H
else
  opt=
fi

$prog --fnhis=histref.dat $opt --fndos=lndosnb0.dat -v > is2nb0.out
$prog --fnhis=histref.dat -H --wham=mdiis --nbases=10 --fndos=lndosnb10.dat -v > is2nb10.out
$prog --fnhis=histref.dat -H --wham=st --fndos=lndosst.dat -v > is2st.out
$prog --fnhis=histref.dat -H --wham=ui --fndos=lndosui.dat -v > is2ui.out

python addlndosref.py

