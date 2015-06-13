#!/usr/bin/env gnuplot

set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced size 10, 4 font "Times, 28"
set output "nsnt.eps"
set multiplot


wl = 0.5
wr = 1 - wl

set rmargin 1.0

ldx = 0.013
ldy = 0.05
set label "(a)" at screen ldx,    1 - ldy font "Times, 32"
set label "(b)" at screen wl+ldx, 1 - ldy font "Times, 32"


set size wl, 1
set origin 0, 0



set xtics 5.0 font "Times,24" offset 0, 0.5
set mxtics 5

set logscale y
set format y "10^{%T}"
set ytics font "Times,24" offset 0.5, 0

color1  = "#224488"
color2  = "#aa0000"
color3  = "#aa00aa"
color4  = "#00aa00"
color5  = "#00aaaa"

set style line 1  lw 2.0 lt 1 lc rgb color1  pt 4   ps 1.8
set style line 2  lw 2.0 lt 2 lc rgb color2  pt 8   ps 2.0
set style line 3  lw 2.0 lt 3 lc rgb color3  pt 10  ps 2.0
set style line 4  lw 2.0 lt 4 lc rgb color4  pt 6   ps 2.0
set style line 5  lw 2.0 lt 5 lc rgb color5  pt 12  ps 2.0

set xlabel "Number of bases, {/Times-Italic M}" offset 0, 0.9
set ylabel "Number of iterations" offset 0.5, 0
set key spacing 1.5
#set key Left reverse spacing 1.5

plot [0:20][:8e3] \
  "../../data/is2/is2wham.dat"    u 1:2 w lp ls 1 t "Ising model", \
  "../../data/1VII/xvgwham.dat"   u 1:2 w lp ls 2 t "Villin headpiece, {/Times-Italic NVT}", \



set size wr, 1
set origin wl, 0

set ylabel "Run time (s)" offset 1.0, 0
#set ytics font "Times,24" offset 0.0, 0.0
#set mytics 4

plot [0:20][:80] \
  "../../data/is2/is2tmwham.dat"    u 1:2 w lp ls 1 t "Ising model", \
  "../../data/1VII/xvgtmwham.dat"   u 1:2 w lp ls 2 t "Villin headpiece, {/Times-Italic NVT}", \



unset multiplot
unset output
set terminal pop
reset
