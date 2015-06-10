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





Notes
======

* MBAR can achieve a lower precision due to the larger number of operations
