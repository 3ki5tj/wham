#!/usr/bin/env gnuplot

set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced size 5, 3.5 font "Times, 28"
set output "nsnt.eps"




set xtics 5.0 font "Times,20" offset 0, 0.5
set mxtics 5
set xlabel "Number of bases, {/Times-Italic M}" offset 0, 0.9

set logscale y
set format y "10^{%T}"
set ytics nomirror font "Times,20" offset 0.5, 0
set ylabel "Number of iterations, {/Times-Italic N}@_{/=16 iter}" offset 0.5, 0
set yrange [:3000]

set logscale y2
set y2tics nomirror font "Times,20" offset -0.5, 0
set y2label "Run time, {/Times-Italic t} (seconds)" textcolor rgb "#404040" offset -3.5, 0
set y2range [:30]
set rmargin 5.0

set key spacing 1.3
#set key Left reverse spacing 1.5

color1  = "#224488"
color2  = "#aa0000"
color3  = "#aa00aa"
color4  = "#00aa00"
color5  = "#00aaaa"

color11 = "#aaaaff"
color12 = "#ffaaaa"
color13 = "#ffaaff"
color14 = "#aaffaa"
color15 = "#aaffff"

set style line 1  lw 2.0 lt 1 lc rgb color1  pt 4   ps 1.8
set style line 2  lw 2.0 lt 2 lc rgb color2  pt 8   ps 2.0
set style line 3  lw 2.0 lt 3 lc rgb color3  pt 10  ps 2.0
set style line 4  lw 2.0 lt 4 lc rgb color4  pt 6   ps 2.0
set style line 5  lw 2.0 lt 5 lc rgb color5  pt 12  ps 2.0

set style line 11 lw 2.0 lt 1 lc rgb color11 pt 5   ps 1.8
set style line 12 lw 2.0 lt 2 lc rgb color12 pt 9   ps 2.0
set style line 13 lw 2.0 lt 3 lc rgb color13 pt 11  ps 2.0
set style line 14 lw 2.0 lt 4 lc rgb color14 pt 7   ps 2.0

plot [0:20][] \
  "../../data/is2/is2tmwham.dat"  u 1:(($1 > -1) ? ($2) : 1/0) axes x1y2 w lp ls 11 t "Ising model, WHAM, {/Times-Italic t}", \
  "../../data/1VII/xvgtmwham.dat" u 1:(($1 > -1) ? ($2) : 1/0) axes x1y2 w lp ls 12 t "Villin headpiece, WHAM, {/Times-Italic t}", \
  "../../data/is2/is2wham.dat"    u 1:(($1 > -1) ? ($2) : 1/0) w lp ls 1 t "Ising model, WHAM, {/Times-Italic N}@_{/=16 iter}", \
  "../../data/1VII/xvgwham.dat"   u 1:(($1 > -1) ? ($2) : 1/0) w lp ls 2 t "Villin headpiece, WHAM, {/Times-Italic N}@_{/=16 iter}", \
  "../../data/1VII/xvgmbar.dat"   u 1:(($1 > -1) ? ($2) : 1/0) w lp ls 3 t "Villin headpiece, MBAR, {/Times-Italic N}@_{/=16 iter}", \



unset output
set terminal pop
reset
