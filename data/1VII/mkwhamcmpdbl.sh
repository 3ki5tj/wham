#!/bin/sh

# make the data files for a comparison
# WHAM, ST-WHAM, MBAR

pdir=../../prog/gmx
prog=$pdir/xvgwham_ldbl
prog2=$pdir/xvgmbar_ldbl

make -C $pdir

#./syncin
#./T_scan.sh

$prog   --wham=mdiis --tol=1e-10 --de=0.1      -v edbl.ls > wham_de0.1dbl.out
$prog   --wham=st                --de=0.1   -H -v edbl.ls > stwham_de0.1dbl.out
$prog   --wham=ui                --de=0.1   -H -v edbl.ls > uiwham_de0.1dbl.out
$prog   --wham=mdiis --tol=1e-10 --de=1        -v edbl.ls --fnhis=histdblde1.dat    > wham_de1dbl.out
$prog   --wham=st                --de=1     -H -v edbl.ls --fnhis=histdblde1.dat    > stwham_de1dbl.out
$prog   --wham=mdiis --tol=1e-10 --de=10       -v edbl.ls --fnhis=histdblde10.dat   > wham_de10dbl.out
$prog   --wham=st                --de=10    -H -v edbl.ls --fnhis=histdblde10.dat   > stwham_de10dbl.out
$prog   --wham=mdiis --tol=1e-10 --de=100      -v edbl.ls --fnhis=histdblde100.dat  > wham_de100dbl.out
$prog   --wham=st                --de=100   -H -v edbl.ls --fnhis=histdblde100.dat  > stwham_de100dbl.out
$prog   --wham=mdiis --tol=1e-10 --de=1000     -v edbl.ls --fnhis=histdblde1000.dat > wham_de1000dbl.out
$prog   --wham=st                --de=1000  -H -v edbl.ls --fnhis=histdblde1000.dat > stwham_de1000dbl.out
$prog2  --mbar=mdiis --tol=1e-10 --itmax=1000  -v edbl.ls > mbardbl.out
$prog2  --est        --tol=1e-10 --itmax=1000  -v edbl.ls > estdbl.out

python addlnzref.py --ref=mbardbl.out *dbl.out

