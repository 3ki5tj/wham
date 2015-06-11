#!/usr/bin/env gnuplot

set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced size 5, 3.5 font "Times, 24"
set output "trace.eps"


set logscale

#set xtics 5.0 font "Times,20" offset 0, 0.5
#set mxtics 5

set format x "10^{%T}"
set xtics font "Times,20" offset 0, 0.0

set format y "10^{%T}"
set ytics font "Times,20" offset 0.5, 0

color1  = "#224488"
color2  = "#aa0000"
color3  = "#aa00aa"
color4  = "#00aa00"
color5  = "#00aaaa"

set style line 1  lw 2.0 lt 1 lc rgb color1  pt 4   ps 0.6
set style line 2  lw 2.0 lt 2 lc rgb color2  pt 8   ps 1.0
set style line 3  lw 2.0 lt 3 lc rgb color3  pt 10  ps 1.0
set style line 4  lw 2.0 lt 4 lc rgb color4  pt 6   ps 1.0
set style line 5  lw 2.0 lt 5 lc rgb color5  pt 12  ps 1.0

set xlabel "Number of iterations" offset 0, 0.0
set ylabel "Error, max\\{{/Times-Italic R_i}\\}" offset 0.5, 0
#set key spacing 1.5
set key left bottom Left reverse spacing 1.5

plot [:2e3][1e-8:] \
  "../../data/is2c/is2kth_nb0wham.dat"   u 1:2 w lp ls 1 t "Direct", \
  "../../data/is2c/is2kth_nb5wham.dat"   u 1:2 w lp ls 2 t "DIIS, {/Times-Italic M} = 5", \
  "../../data/is2c/is2kth_nb10wham.dat"  u 1:2 w lp ls 3 t "DIIS, {/Times-Italic M} = 10", \
  "../../data/is2c/is2kth_nb15wham.dat"  u 1:2 w lp ls 4 t "DIIS, {/Times-Italic M} = 15", \



unset output
set terminal pop
reset
