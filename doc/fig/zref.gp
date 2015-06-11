#!/usr/bin/env gnuplot

set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced size 5, 3.5 font "Times, 24"
set output "zref.eps"
set multiplot


set xtics 0.5 font "Times,20" offset 0, 0.5
set mxtics 5

set ytics 1000 font "Times,20" offset 0, 0
set mytics 5


color1  = "#224488"
color2  = "#aa0000"
color3  = "#aa00aa"
color4  = "#00aa00"
color9  = "#000000"

set style line 1  lw 2.0 lt 1 lc rgb color1  pt 4   ps 0.8
set style line 2  lw 2.0 lt 2 lc rgb color2  pt 8   ps 1.0
set style line 3  lw 2.0 lt 3 lc rgb color3  pt 10  ps 1.0
set style line 4  lw 2.0 lt 4 lc rgb color4  pt 6   ps 1.0

set style line 9  lw 2.0 lt 1 lc rgb color9  pt 1   ps 1.0

set xlabel "{/Times-Italic T_i}" offset 0, 0.9
set ylabel "{/Times-Italic f_i} = -log {/Times-Italic Z_i}" offset 2.0, 0
set key left Left reverse spacing 1.5


plot [1.5:3.1][:2400] \
  "../../data/is2/is2nb0.out"   u (1/$2):(-$3) w p ls 1 t "Direct WHAM", \
  "../../data/is2/is2nb10.out"  u (1/$2):(-$3) w p ls 2 t "DIIS WHAM, {/Times-Italic M} = 10", \
  ""                            u (1/$2):(-$4) w l ls 9 t "Reference", \



set size 0.55, 0.5
set origin 0.4, 0.15

unset xlabel
set ylabel "{/Times-Italic f_i} - {/Times-Italic f_i}^{/Times (ref)}" offset 3.0, 1.0
set ytics 0.1 font "Times, 16" offset 0.5, 0
set xtics font "Times, 16"
set format x ""


plot [1.5:3.1][-0.18:] \
  "../../data/is2/is2nb0.out"   u (1/$2):(-$3+$4) w p ls 1 notitle, \
  "../../data/is2/is2nb10.out"  u (1/$2):(-$3+$4) w p ls 2 notitle, \


unset multiplot
unset output
set terminal pop
reset
