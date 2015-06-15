#!/bin/sh

# make the data files for a comparison
# WHAM, ST-WHAM, MBAR

pdir=../../prog/gmx
prog=$pdir/xvgwham
prog2=$pdir/xvgmbar

make -C $pdir

#./syncin
#./T_scan.sh

$prog   --wham=mdiis --de=0.1      -v e.ls > wham_de0.1.out
$prog   --wham=st    --de=0.1   -H -v e.ls > stwham_de0.1.out
$prog   --wham=mdiis --de=1        -v e.ls > wham_de1.out
$prog   --wham=st    --de=1     -H -v e.ls > stwham_de1.out
$prog   --wham=mdiis --de=3        -v e.ls > wham_de3.out
$prog   --wham=st    --de=3     -H -v e.ls > stwham_de3.out
$prog   --wham=mdiis --de=10       -v e.ls > wham_de10.out
$prog   --wham=st    --de=10    -H -v e.ls > stwham_de10.out
$prog   --wham=mdiis --de=30       -v e.ls > wham_de30.out
$prog   --wham=st    --de=30    -H -v e.ls > stwham_de30.out
$prog   --wham=mdiis --de=100      -v e.ls > wham_de100.out
$prog   --wham=st    --de=100   -H -v e.ls > stwham_de100.out
$prog   --wham=mdiis --de=300      -v e.ls > wham_de300.out
$prog   --wham=st    --de=300   -H -v e.ls > stwham_de300.out
$prog   --wham=mdiis --de=1000     -v e.ls > wham_de1000.out
$prog   --wham=st    --de=1000  -H -v e.ls > stwham_de1000.out
$prog2  --mbar=mdiis --tol=2e-8    -v e.ls > mbar.out

python addlnzref.py
