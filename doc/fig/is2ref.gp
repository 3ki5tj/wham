#!/usr/bin/env gnuplot



set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced size 10, 4.5 font "Times, 40"
set output "is2ref.eps"
set multiplot


wl = 0.5
wr = 1 - wl
ht = 0.85


n = 64*64

set bmargin 2.8
set rmargin 1.0


ldx = 0.013
ldy = 0.06
set label "(a)" at screen ldx,    1 - ldy font "Times, 48"
set label "(b)" at screen wl+ldx, 1 - ldy font "Times, 48"



set size wl, ht
set origin 0, 1 - ht



color1  = "#aa0000"
color2  = "#224488"
color3  = "#006000"
color4  = "#aa00aa"
color9  = "#808080"

set style line 1  lw 2.0 lt 1 lc rgb color1  pt 10  ps 1.6
set style line 2  lw 2.0 lt 1 lc rgb color2  pt 8   ps 1.6
set style line 3  lw 2.0 lt 1 lc rgb color3  pt 4   ps 1.2
set style line 4  lw 2.0 lt 1 lc rgb color4  pt 6   ps 1.2

set style line 9  lw 0.5 lt 4 lc rgb color9  pt 1   ps 1.0

set xlabel "{/Times-Italic T_i}" offset 0, 0.5
set xtics 0.5 offset 0, 0
set mxtics 5

set ylabel "{/Symbol-Oblique e} ({/Times-Italic f_i})" offset 1.5, 0.0
set ytics 0.1 offset 0.5, 0
set mytics 5

set key at screen 0.0, 0.03 left bottom Left reverse width -6 maxrows 1 samplen 3



plot [1.5:3.1][-0.14:0.02] \
  0 w l ls 9 notitle, \
  "../../data/is2/is2nb0.out"   u (1/$2):(-$3+$4) w p ls 1 notitle, \
  "../../data/is2/is2nb10.out"  u (1/$2):(-$3+$4) w p ls 2 notitle, \
  "../../data/is2/is2st.out"    u (1/$2):(-$3+$4) w p ls 3 notitle, \
  "../../data/is2/is2ui.out"    u (1/$2):(-$3+$4) w p ls 4 notitle, \
  -10 w lp ls 1 ps 3.0 t "Direct WHAM", \
  -10 w lp ls 2 ps 3.0 t "DIIS WHAM, {/Times-Italic M} = 10"




set size wr, ht
set origin wl, 1 - ht

set xlabel "{/Times-Italic E} / {/Times-Italic N}" offset 2, 0.5
set xtics 0.5 offset 0, 0
set mxtics 5
set format x "%g"

set ylabel "{/Symbol-Oblique e} [log {/Times-Italic g}({/Times-Italic E})]" offset 1.5, 0.0
set ytics 0.1 offset 0.5, 0
set mytics 5

set key at screen 0.61, 0.03 width +2

plot [:-0.7][-0.14:0.1] \
  0 w l ls 9 notitle, \
  "../../data/is2/lndosnb0.dat"   u ($1/n):($2-$3) w l ls 1 notitle, \
  "../../data/is2/lndosnb10.dat"  u ($1/n):($2-$3) w l ls 2 notitle, \
  "../../data/is2/lndosst.dat"    u ($1/n):($2-$3) w l ls 3 notitle, \
  "../../data/is2/lndosui.dat"    u ($1/n):($2-$3) w l ls 4 notitle, \
  "../../data/is2/lndosnb0.dat"   u ($1/n):($2-$3) every 50 w p ls 1 ps 2.0 notitle, \
  "../../data/is2/lndosnb10.dat"  u ($1/n):($2-$3) every 50 w p ls 2 ps 2.0 notitle, \
  "../../data/is2/lndosst.dat"    u ($1/n):($2-$3) every 50 w p ls 3 ps 1.6 notitle, \
  "../../data/is2/lndosui.dat"    u ($1/n):($2-$3) every 50 w p ls 4 ps 1.6 notitle, \
  -10 w lp ls 3 ps 2.3 t "ST-WHAM", \
  -10 w lp ls 4 ps 2.3 t "UIM"



unset multiplot
unset output
set terminal pop
reset
