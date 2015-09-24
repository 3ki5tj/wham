#!/usr/bin/env gnuplot

set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced size 10, 7 font "Times, 48"
set output "is2approx.eps"


set multiplot


n = 1024

wl = 0.54
wr = 1 - wl
ht = 0.57
ht2 = 0.35

set size wl, ht
set origin 0, 1 - ht



set xlabel "{/Times-Italic T}" offset -1.0, 1.0
set xtics 0.5 offset 0, 0.3
set mxtics 5
set xrange [1.5:3.2]

set logscale y
set format y "10^{%T}"
set ylabel "&{i}_{/*2.5 |} {/Times-Itaic f_i} - {/Times-Italic f_i}^{ref}&{i}_{/*2.5 |}" offset 0, 0
set ytics offset 0.5, 0
set mytics 10
set yrange [2e-6:40]

set tmargin 1
set lmargin 7
set rmargin 0

color1  = "#882222"
color2  = "#224488"
color3  = "#8822aa"
color4  = "#006000"
color5  = "#008888"

color11 = "#ffaaaa"
color12 = "#aaaaff"
color13 = "#ddaaff"
color14 = "#80cc80"



set style line 1  lw 1.0 lt 4 lc rgb color1  pt 5  ps 2.3
set style line 2  lw 1.0 lt 4 lc rgb color2  pt 9  ps 3.0
set style line 3  lw 1.0 lt 4 lc rgb color3  pt 11 ps 3.0
set style line 4  lw 1.0 lt 4 lc rgb color4  pt 7  ps 2.5
set style line 5  lw 1.0 lt 4 lc rgb color5  pt 13 ps 3.0

set key at screen 0, ht2 + 0.02 left bottom Left reverse invert samplen 4 width 0.5 maxrows 1

set label 1 "(a) {/Symbol-Oblique D}{/Times-Italic T} = 0.05" at 1.55, 8.0

plot \
  "../../data/is2/is2approx0.05.out" u 1:($1 > 1.501 ? abs($2-$3) : 1/0) w lp ls 1 notitle, \
  "../../data/is2/is2approx0.05.out" u 1:($1 > 1.501 ? abs($2-$4) : 1/0) w lp ls 2 notitle, \
  "../../data/is2/is2approx0.05.out" u 1:($1 > 1.501 ? abs($2-$5) : 1/0) w lp ls 3 notitle, \
  "../../data/is2/is2approx0.05.out" u 1:($1 > 1.501 ? abs($2-$9) : 1/0) w lp ls 5 notitle, \
  "../../data/is2/is2approx0.05.out" u 1:($1 > 1.501 ? abs($2-$7) : 1/0) w lp ls 4 notitle, \
  1e10 w lp ls 1 ps 3.2 t "TPA", \
  1e10 w lp ls 2 ps 4.0 t "EM", \
  1e10 w lp ls 3 ps 4.0 t "VEM", \
  1e10 w lp ls 5 ps 4.0 t "TG", \
  1e10 w lp ls 4 ps 3.2 t "UIM"



set size wr, ht
set origin wl, 1 - ht


unset ylabel
set format y ""
set lmargin 0
set rmargin 1

set label 1 "(b) {/Symbol-Oblique D}{/Times-Italic T} = 0.2"

plot \
  "../../data/is2/is2approx0.2.out" u 1:($1 > 1.501 ? abs($2-$3) : 1/0) w lp ls 1 ps 3.2 notitle, \
  "../../data/is2/is2approx0.2.out" u 1:($1 > 1.501 ? abs($2-$4) : 1/0) w lp ls 2 ps 4.0 notitle, \
  "../../data/is2/is2approx0.2.out" u 1:($1 > 1.501 ? abs($2-$5) : 1/0) w lp ls 3 ps 4.0 notitle, \
  "../../data/is2/is2approx0.2.out" u 1:($1 > 1.501 ? abs($2-$9) : 1/0) w lp ls 5 ps 4.0 notitle, \
  "../../data/is2/is2approx0.2.out" u 1:($1 > 1.501 ? abs($2-$7) : 1/0) w lp ls 4 ps 3.2 notitle, \


reset

set size 1.0, ht2
set origin 0, 0

set xtics 0.5 offset 0, 0.5
set mxtics 5

unset ytics
#set ytics 0.05
#set mytics 5

set tmargin 0.5
set lmargin 7
set rmargin 1.0

set label 1 at -1.5, 0.08 "(c) {/Symbol-Oblique D}{/Times-Italic T} = 0.2"

set xlabel "{/Times-Italic E}/{/Times-Italic N}" offset 0.0, 1.0

set ylabel "Energy\nhistograms,\n{/Times-Itaic H}({/Times-Italic E})" offset 0, 0

plot [:-0.6] \
  "../../data/is2/histref32x32.dat" u ($1/n):($2) w l lw 3 notitle



unset multiplot

unset output
set terminal pop
reset
