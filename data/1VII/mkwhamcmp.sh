#!/bin/sh

# make the data files for a comparison
# WHAM, ST-WHAM, MBAR

pdir=../../prog/gmx
prog=$pdir/xvgwham
prog2=$pdir/xvgmbar

make -C $pdir

#./syncin
#./T_scan.sh

$prog   --wham=mdiis --tol=1e-9 --de=0.1      -v e.ls > wham_de0.1.out
$prog   --wham=st               --de=0.1   -H -v e.ls > stwham_de0.1.out
$prog   --wham=ui               --de=0.1   -H -v e.ls > uiwham_de0.1.out
$prog   --wham=mdiis --tol=1e-9 --de=1        -v e.ls > wham_de1.out
$prog   --wham=st               --de=1     -H -v e.ls > stwham_de1.out
$prog   --wham=mdiis --tol=1e-9 --de=10       -v e.ls > wham_de10.out
$prog   --wham=st               --de=10    -H -v e.ls > stwham_de10.out
$prog   --wham=mdiis --tol=1e-9 --de=100      -v e.ls > wham_de100.out
$prog   --wham=st               --de=100   -H -v e.ls > stwham_de100.out
$prog   --wham=mdiis --tol=1e-9 --de=1000     -v e.ls > wham_de1000.out
$prog   --wham=st               --de=1000  -H -v e.ls > stwham_de1000.out
$prog2  --mbar=mdiis --tol=2.5e-8 --itmax=100 -v e.ls > mbar.out
$prog2  --est                                 -v e.ls > est.out

python addlnzref.py
