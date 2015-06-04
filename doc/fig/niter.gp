#!/usr/bin/env gnuplot

set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced size 5, 3.5 font "Times, 24"
set output "niter.eps"



set xtics 5.0 font "Times,20" offset 0, 0.5
set mxtics 5

set logscale y
set format y "10^{%T}"
set ytics font "Times,20" offset 0.5, 0

color1  = "#224488"
color2  = "#aa0000"
color3  = "#aa00aa"

set style line 1  lw 2.0 lt 1 lc rgb color1  pt 4   ps 1.8
set style line 2  lw 2.0 lt 2 lc rgb color2  pt 8   ps 2.0
set style line 3  lw 2.0 lt 3 lc rgb color3  pt 10  ps 2.0

set xlabel "Number of bases" offset 0, 0.9
set ylabel "Number of iterations" offset 0.5, 0
set key spacing 1.5
#set key Left reverse spacing 1.5

plot [0:20][:] \
  "../../data/is2/is2wham.dat"  u 1:(($1 > 0) ? $2 : 1/0) w p ls 1 t "Ising model", \
  "../../data/lj/ljwham.dat"    u 1:(($1 > 0) ? $2 : 1/0) w p ls 2 t "Lennard-Jones fluid, {/Times-Italic NVT}", \
  "../../data/lj/lj2wham.dat"   u 1:(($1 > 0) ? $2 : 1/0) w p ls 3 t "Lennard-Jones fluid, {/Times-Italic NpT}", \



unset output
set terminal pop
reset
