#!/bin/sh

# make the data files for is2approx.gp

pdir=../../prog/is2
prog=$pdir/is2ref
side=32

make -C $pdir $prog

$prog --dT=0.05  --nT=33  --L=$side -v > is2approx0.05.out
$prog --dT=0.10  --nT=17  --L=$side -v > is2approx0.1.out
$prog --dT=0.20  --nT=9   --L=$side -v > is2approx0.2.out


