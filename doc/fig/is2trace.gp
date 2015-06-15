#!/usr/bin/env gnuplot

set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced size 5, 3.5 font "Times, 24"
set output "is2trace.eps"


set logscale

set format x "10^{%T}"
set xtics font "Times,20" offset 0, 0.0
set mxtics 10

set format y "10^{%T}"
set ytics font "Times,20" offset 0.5, 0
set mytics 10

color1  = "#224488"
color2  = "#aa0000"
color3  = "#aa00aa"
color4  = "#00aa00"
color5  = "#00aaaa"

set style line 1  lw 1.0 lt 1 lc rgb color1  pt 1   ps 0.6
set style line 2  lw 1.0 lt 2 lc rgb color2  pt 8   ps 1.0
set style line 3  lw 1.0 lt 3 lc rgb color3  pt 10  ps 1.0
set style line 4  lw 1.0 lt 4 lc rgb color4  pt 6   ps 1.0
set style line 5  lw 1.0 lt 5 lc rgb color5  pt 12  ps 1.0

set xlabel "Number of iterations" offset 0, 0.0
set ylabel "max{/Time=36 \\{|}{/Times-Italic R_i}{/Time=36 |\\}}" offset 0.5, 0
#set key spacing 1.5
set key left bottom Left reverse spacing 1.5 at 0.8, 2e-8

plot [:3e3][1e-8:] \
  "../../data/is2/is2_nb0wham.dat"   u 1:2 w p ls 1 t "Direct", \
  "../../data/is2/is2_nb5wham.dat"   u 1:2 w p ls 2 t "DIIS, {/Times-Italic M} = 5", \
  "../../data/is2/is2_nb10wham.dat"  u 1:2 w p ls 3 t "DIIS, {/Times-Italic M} = 10", \
  "../../data/is2/is2_nb15wham.dat"  u 1:2 w p ls 4 t "DIIS, {/Times-Italic M} = 15", \



unset output
set terminal pop
reset
