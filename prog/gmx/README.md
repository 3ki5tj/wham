Overview
========

Weighted histogram analysis method (WHAM)
and multistate Bennett's acceptance ratio (MBAR) method
applied on .xvg files (output of GROMACS)



Files
=====

File        | Description
------------|-----------------------------------
xvgact.c    | Autocorrelation time
xvgwham.c   | NVT WHAM driver
xvgwham2.c  | NPT WHAM driver
xvgmbar.c   | NVT MBAR driver
xvgmbar2.c  | NPT MBAR driver
mbar.h      | NVT MBAR core file
mbar2.h     | NPT MBAR core file
lsutil.h    | utility to handle a list of NVT file names
ls2util.h   | utility to handle a list of NPT file names



xdouble versions
=================

WHAM and MBAR programs can be compiled with "-DLDBL" for long double precision.



Notes
======

* MBAR can achieve a lower precision due to the larger number of operations



Usage
======


### xvgact

```
./xvgact ../data/1VII/T300/e.xvg --fnac=ac.dat
```

### xvgwham


Direct WHAM, simplest example.
```
./xvgwham e.ls
```

MDIIS-WHAM, be verbose
```
./xvgwham e.ls --wham=mdiis -v
```
The options `--wham=` can also be `st` and `ui`.


MDIIS-WHAM, save a histogram and a file for autocorrelation times,
and then do bootscrapping three times
```
./xvgwham e.ls --fnhis=hist.dat --wham=mdiis --fnact=act.dat
./xvgwham e.ls --fnhis=hist.dat --wham=mdiis --fnact=act.dat --bootstrap -H
./xvgwham e.ls --fnhis=hist.dat --wham=mdiis --fnact=act.dat --bootstrap -H
./xvgwham e.ls --fnhis=hist.dat --wham=mdiis --fnact=act.dat --bootstrap -H
```



