#!/bin/sh

# make the data files for is2approx.gp

pdir=../../prog/is2
prg=is2ref
prghis=is2mkhis
prog=$pdir/$prg
side=32

make -C $pdir $prg $prghis

$prog --dT=0.02  --nT=81  --L=$side -v > is2approx0.02.out
$prog --dT=0.05  --nT=33  --L=$side -v > is2approx0.05.out
$prog --dT=0.10  --nT=17  --L=$side -v > is2approx0.1.out
$prog --dT=0.20  --nT=9   --L=$side -v > is2approx0.2.out

$pdir/$prghis
