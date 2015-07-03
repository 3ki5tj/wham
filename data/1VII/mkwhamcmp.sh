#!/bin/sh

# make the data files for a comparison
# WHAM, ST-WHAM, MBAR

pdir=../../prog/gmx
prog=$pdir/xvgwham_ldbl
prog2=$pdir/xvgmbar_ldbl

make -C $pdir

#./syncin
#./T_scan.sh

$prog   --wham=mdiis --tol=1e-10 --de=0.1      -v e.ls > wham_de0.1.out
$prog   --wham=st                --de=0.1   -H -v e.ls > stwham_de0.1.out
$prog   --wham=ui                --de=0.1   -H -v e.ls > uiwham_de0.1.out
$prog   --wham=mdiis --tol=1e-10 --de=1        -v e.ls --fnhis=histde1.dat    > wham_de1.out
$prog   --wham=st                --de=1     -H -v e.ls --fnhis=histde1.dat    > stwham_de1.out
$prog   --wham=mdiis --tol=1e-10 --de=10       -v e.ls --fnhis=histde10.dat   > wham_de10.out
$prog   --wham=st                --de=10    -H -v e.ls --fnhis=histde10.dat   > stwham_de10.out
$prog   --wham=mdiis --tol=1e-10 --de=100      -v e.ls --fnhis=histde100.dat  > wham_de100.out
$prog   --wham=st                --de=100   -H -v e.ls --fnhis=histde100.dat  > stwham_de100.out
$prog   --wham=mdiis --tol=1e-10 --de=1000     -v e.ls --fnhis=histde1000.dat > wham_de1000.out
$prog   --wham=st                --de=1000  -H -v e.ls --fnhis=histde1000.dat > stwham_de1000.out
$prog2  --mbar=mdiis --tol=1e-10 --itmax=1000  -v e.ls > mbar.out
$prog2  --est        --tol=1e-10 --itmax=1000  -v e.ls > est.out

python addlnzref.py
