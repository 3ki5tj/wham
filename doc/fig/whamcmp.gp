#!/usr/bin/env gnuplot



set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced size 5, 4.5 font "Times, 24"
set output "whamcmp.eps"

set multiplot


ht1 = 0.72
ht2 = 1 - ht1

set label 1 "(a)" at screen 0.01, 1 - 0.02   font "Times, 32"
set label 2 "(b)" at screen 0.01, ht2 + 0.0  font "Times, 32"


set size 1.0, ht1
set origin 0.0, ht2


BOLTZ = 1.380658*6.0221367/1000

set xlabel "{/Times-Italic T_i} (K)" offset 1, 1
set xtics 50 font "Times,20" offset 0, 0.5
set mxtics 5

set logscale y
set format y "10^{%T}"
set ytics font "Times,20" offset 0.5, 0
set mytics 20
set ylabel "&{i}_{/*2.4|} {/Times-Italic f_i} - {/Times-Italic f_i}^{/*0.7 (MBAR)}&{i}_{/*2.4|}" offset 1.0, 1.0

set key at 433, 2e-5 Left reverse font "Times,18" width -6 vertical maxrows 3 samplen 2.5 spacing 1.0

color1  = "#224488"
color2  = "#aa2222"
color3  = "#804000"
color4  = "#aa00aa"
color5  = "#208020"
color6  = "#6688ff"
color7  = "#228888"

color9  = "#808080"

set style line 1  lw 2.0 lt 4 lc rgb color1  pt 8   ps 2.0
set style line 2  lw 2.0 lt 4 lc rgb color1  pt 5   ps 1.6
set style line 3  lw 2.0 lt 4 lc rgb color1  pt 13  ps 2.0
set style line 4  lw 2.0 lt 4 lc rgb color6  pt 7   ps 1.8

set style line 11 lw 2.0 lt 4 lc rgb color2  pt 9   ps 2.0
set style line 12 lw 2.0 lt 4 lc rgb color2  pt 11  ps 1.6
set style line 13 lw 2.0 lt 4 lc rgb color4  pt 15  ps 1.8
set style line 14 lw 2.0 lt 4 lc rgb color2  pt 13  ps 2.0

set style line 21 lw 3.0 lt 4 lc rgb color3  pt 1   ps 1.6
set style line 22 lw 3.0 lt 4 lc rgb color5  pt 2   ps 1.2
set style line 23 lw 3.0 lt 4 lc rgb color7  pt 3   ps 1.2
set style line 24 lw 3.0 lt 4 lc rgb color7  pt 15  ps 1.6

set style line 9  lw 2.0 lt 1 lc rgb color9  pt 10  ps 1.5

plot [300:420][1e-6:1] \
  "../../data/1VII/wham_de1.out"      u (1/BOLTZ/$2):(($1 > 0) ? abs($7-$3)  : 1/0) w lp ls 2   t "WHAM, {/Times-Italic h} = 1", \
  "../../data/1VII/wham_de100.out"    u (1/BOLTZ/$2):(($1 > 0) ? abs($7-$3)  : 1/0) w lp ls 3   t "WHAM, {/Times-Italic h} = 100", \
  "../../data/1VII/est.out"           u (1/BOLTZ/$2):(($1 > 0) ? abs($12-$8) : 1/0) w lp ls 4   t "BAR, ({/Times-Italic i}, {/Times-Italic i} + 1)", \
  "../../data/1VII/stwham_de1.out"    u (1/BOLTZ/$2):(($1 > 0) ? abs($7-$3)  : 1/0) w lp ls 11  t "ST-WHAM, {/Times-Italic h} = 1", \
  "../../data/1VII/stwham_de100.out"  u (1/BOLTZ/$2):(($1 > 0) ? abs($7-$3)  : 1/0) w lp ls 12  t "ST-WHAM, {/Times-Italic h} = 100", \
  "../../data/1VII/uiwham_de0.1.out"  u (1/BOLTZ/$2):(($1 > 0) ? abs($7-$3)  : 1/0) w lp ls 13  t "UIM", \
  "../../data/1VII/est.out"           u (1/BOLTZ/$2):(($1 > 0) ? abs($12-$3) : 1/0) w lp ls 21  t "Eq. (10)", \
  "../../data/1VII/est.out"           u (1/BOLTZ/$2):(($1 > 0) ? abs($12-$4) : 1/0) w lp ls 22  t "Eq. (11)", \


# "../../data/1VII/est.out"           u (1/BOLTZ/$2):(($1 > 0) ? abs($12-$5) : 1/0) w lp ls 23  t "Eq. (C3)", \
#  "../../data/1VII/xvgbswham.dat"     u (1/BOLTZ/$2):(($1 > 0) ? abs($4)     : 1/0) w l  ls 9   t "Random error of WHAM",


reset

set size 1.0, ht2
set origin 0, 0

set tmargin 0
set lmargin 6

set xlabel "{/Times-Italic E}" offset -1, 0.5
set xtics 5000
set mxtics 5

set ylabel "Energy\nhistograms,\n{/Times-Italic H}({/Times-Italic E})" offset 0, -0.5
unset ytics

plot [-82000:-70000][] "../../data/1VII/histde1.dat" u 1:2 w l notitle



unset multiplot
unset output
set terminal pop
reset
