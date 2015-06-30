#!/usr/bin/env gnuplot

set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced size 10, 4 font "Times, 24"
set output "is2approx.eps"


set multiplot


wl = 0.54
wr = 1 - wl



set size wl, 1


set xlabel "{/Times-Italic T}" font "Times,32" offset 0, 0.5
set xtics 0.5 font "Times,28" offset 0, 0.3
set mxtics 5

set logscale y
set format y "10^{%T}"
set ylabel "&{i}_{/*2.5 |} {/Times-Itaic f_i} - {/Times-Italic f_i}^{ref}&{i}_{/*2.5 |}" font "Times,32" offset -1, 0
set ytics font "Times,28" offset 0.5, 0
set mytics 10

set rmargin 0.

color1  = "#aa2222"
color2  = "#224488"
color3  = "#8822aa"
color4  = "#008000"

color11 = "#ffaaaa"
color12 = "#aaaaff"
color13 = "#ddaaff"
color14 = "#80cc80"



set style line 1  lw 1.0 lt 1 lc rgb color1  pt 4   ps 1.6
set style line 2  lw 1.0 lt 2 lc rgb color2  pt 8   ps 2.0
set style line 3  lw 1.0 lt 3 lc rgb color3  pt 10  ps 2.0
set style line 4  lw 1.0 lt 4 lc rgb color4  pt 6   ps 2.0

set style line 11 lw 1.0 lt 1 lc rgb color11 pt 5   ps 1.6
set style line 12 lw 1.0 lt 2 lc rgb color12 pt 9   ps 2.0
set style line 13 lw 1.0 lt 3 lc rgb color13 pt 11  ps 2.0
set style line 14 lw 1.0 lt 4 lc rgb color14 pt 7   ps 2.0

set style line 21 lw 1.0 lt 1 lc rgb color1  pt 5   ps 1.6
set style line 22 lw 1.0 lt 2 lc rgb color2  pt 9   ps 2.0
set style line 23 lw 1.0 lt 3 lc rgb color3  pt 11  ps 2.0
set style line 24 lw 1.0 lt 4 lc rgb color4  pt 7   ps 2.0

set key at 3.1, 4e-6 Left reverse invert maxrows 4 samplen 1.8 width -8.0 font "Times,21" spacing 1.3

set label 1 "(a)" at 1.55, 1.2 font "Times,48"

plot [1.5:3.1][1e-7:3] \
  "../../data/is2/is2approx0.2.out"  u 1:($1 > 1.501 ? abs($2-$3) : 1/0) w lp ls 21 t "Eq. (C1), {/Symbol-Oblique D}{/Times-Italic T} = 0.2", \
  "../../data/is2/is2approx0.1.out"  u 1:($1 > 1.501 ? abs($2-$3) : 1/0) w lp ls 11 t "Eq. (C1), {/Symbol-Oblique D}{/Times-Italic T} = 0.1", \
  "../../data/is2/is2approx0.05.out" u 1:($1 > 1.501 ? abs($2-$3) : 1/0) w lp ls  1 t "Eq. (C1), {/Symbol-Oblique D}{/Times-Italic T} = 0.05", \
  "../../data/is2/is2approx0.2.out"  u 1:($1 > 1.501 ? abs($2-$4) : 1/0) w lp ls 22 t "Eq. (C2), {/Symbol-Oblique D}{/Times-Italic T} = 0.2", \
  "../../data/is2/is2approx0.1.out"  u 1:($1 > 1.501 ? abs($2-$4) : 1/0) w lp ls 12 t "Eq. (C2), {/Symbol-Oblique D}{/Times-Italic T} = 0.1", \
  "../../data/is2/is2approx0.05.out" u 1:($1 > 1.501 ? abs($2-$4) : 1/0) w lp ls  2 t "Eq. (C2), {/Symbol-Oblique D}{/Times-Italic T} = 0.05", \
  "../../data/is2/is2approx0.2.out"  u 1:($1 > 1.501 ? abs($2-$7) : 1/0) w lp ls 24 t "UIM, {/Symbol-Oblique D}{/Times-Italic T} = 0.2", \
  "../../data/is2/is2approx0.1.out"  u 1:($1 > 1.501 ? abs($2-$7) : 1/0) w lp ls 14 t "UIM, {/Symbol-Oblique D}{/Times-Italic T} = 0.1", \
  "../../data/is2/is2approx0.05.out" u 1:($1 > 1.501 ? abs($2-$7) : 1/0) w lp ls  4 t "UIM, {/Symbol-Oblique D}{/Times-Italic T} = 0.05", \


set origin wl, 0
set size wr, 1


unset ylabel
set format y ""
set lmargin 0
set rmargin 1

set key at 2.9, 4e-6

set label 1 "(b)"

plot [1.5:3.1][1e-7:3] \
  "../../data/is2/is2approx0.2.out"  u 1:($1 > 1.501 ? abs($2-$4) : 1/0) w lp ls 22 t "Eq. (C2), {/Symbol-Oblique D}{/Times-Italic T} = 0.2", \
  "../../data/is2/is2approx0.1.out"  u 1:($1 > 1.501 ? abs($2-$4) : 1/0) w lp ls 12 t "Eq. (C2), {/Symbol-Oblique D}{/Times-Italic T} = 0.1", \
  "../../data/is2/is2approx0.05.out" u 1:($1 > 1.501 ? abs($2-$4) : 1/0) w lp ls  2 t "Eq. (C2), {/Symbol-Oblique D}{/Times-Italic T} = 0.05", \
  "../../data/is2/is2approx0.2.out"  u 1:($1 > 1.501 ? abs($2-$5) : 1/0) w lp ls 23 t "Eq. (C3), {/Symbol-Oblique D}{/Times-Italic T} = 0.2", \
  "../../data/is2/is2approx0.1.out"  u 1:($1 > 1.501 ? abs($2-$5) : 1/0) w lp ls 13 t "Eq. (C3), {/Symbol-Oblique D}{/Times-Italic T} = 0.1", \
  "../../data/is2/is2approx0.05.out" u 1:($1 > 1.501 ? abs($2-$5) : 1/0) w lp ls  3 t "Eq. (C3), {/Symbol-Oblique D}{/Times-Italic T} = 0.05", \


unset multiplot

unset output
set terminal pop
reset
