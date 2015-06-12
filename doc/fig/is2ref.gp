#!/usr/bin/env gnuplot



set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced size 10, 4 font "Times, 28"
set output "is2ref.eps"
set multiplot


wl = 0.5
wr = 1 - wl

n = 64*64

set bmargin 3
set rmargin 0.5


ldx = 0.013
ldy = 0.05
set label "(a)" at screen ldx,    1 - ldy font "Times, 32"
set label "(b)" at screen wl+ldx, 1 - ldy font "Times, 32"


set size wl, 1
set origin 0, 0



set xtics 0.5 font "Times,24" offset 0, 0.0
set mxtics 5

set ytics 0.2 font "Times,24" offset 0, 0
set mytics 4


color1  = "#224488"
color2  = "#aa0000"
color3  = "#00aa00"
color4  = "#aa00aa"
color9  = "#000000"

set style line 1  lw 1.0 lt 1 lc rgb color1  pt 4   ps 0.8
set style line 2  lw 1.0 lt 2 lc rgb color2  pt 8   ps 1.0
set style line 3  lw 1.0 lt 3 lc rgb color3  pt 6   ps 0.8
set style line 4  lw 1.0 lt 4 lc rgb color4  pt 10  ps 1.0

set style line 9  lw 1.0 lt 1 lc rgb color9  pt 1   ps 1.0

set xlabel "{/Times-Italic T_i}" offset 0, 0.0
set ylabel "{/Times-Italic f_i} / {/Times-Italic N} = -[ log {/Times-Italic Z_i} ] / {/Times-Italic N}" offset 2.0, 0
set key left Left reverse spacing 1.5


plot [1.5:3.1][:] \
  "../../data/is2/is2nb0.out"   u (1/$2):(-$3/n) w p ls 1 t "Direct WHAM", \
  "../../data/is2/is2nb10.out"  u (1/$2):(-$3/n) w p ls 2 t "DIIS WHAM, {/Times-Italic M} = 10", \
  "../../data/is2/is2st.out"    u (1/$2):(-$3/n) w p ls 3 t "ST-WHAM", \
  "../../data/is2/is2ref.dat"   u ($1):(-$2/n)   w l ls 9 t "Reference", \




# inset for f
set size 0.55*wl, 0.6
set origin 0.4*wl, 0.05

unset xlabel
set xtics font "Times, 16"
set format x ""

set ylabel "{/Symbol-Oblique e} ({/Times-Italic f_i})" offset 3.0, 0.0
set ytics 0.1 font "Times, 16" offset 0.5, 0
set mytics 5


plot [1.5:3.1][:] \
  "../../data/is2/is2nb0.out"   u (1/$2):(-$3+$4) w p ls 1 notitle, \
  "../../data/is2/is2nb10.out"  u (1/$2):(-$3+$4) w p ls 2 notitle, \
  "../../data/is2/is2st.out"    u (1/$2):(-$3+$4) w p ls 3 notitle, \




set size wr, 1
set origin wl, 0

set xlabel "{/Times-Italic E} / {/Times-Italic N}" offset 0, 0.0
set xtics 0.5 font "Times,24" offset 0.0, 0.0
set mxtics 5
set format x "%g"

set ylabel "[ log {/Times-Italic g}({/Times-Italic E}) ] / {/Times-Italic N}" offset 2.0, 0.0
set ytics 0.2 font "Times,24" offset 0.0, 0.0
set mytics 4

plot [:-0.7][:] \
  "../../data/is2/lndosnb0.dat"   u ($1/n):($2/n) every 20 w p ls 1 t "Direct WHAM", \
  "../../data/is2/lndosnb10.dat"  u ($1/n):($2/n) every 20 w p ls 2 t "DIIS WHAM", \
  "../../data/is2/lndosst.dat"    u ($1/n):($2/n) every 20 w p ls 3 t "ST WHAM", \
  ""                              u ($1/n):($3/n) w l ls 9 t "Reference", \



# inset for ln g(E)
set size 0.55*wr, 0.5
set origin wl + 0.4*wr, 0.05

unset xlabel
set xtics font "Times, 16"
set format x ""

set ylabel "{/Symbol-Oblique e} [ log {/Times-Italic g}({/Times-Italic E}) ]" offset 3.0, 0.0
set ytics 0.1 font "Times, 16" offset 0.5, 0
set mytics 5


plot [:-0.7][-0.16:0.36] \
  "../../data/is2/lndosnb0.dat"   u ($1/n):($2-$3) every 20 w p ls 1 notitle, \
  "../../data/is2/lndosnb10.dat"  u ($1/n):($2-$3) every 20 w p ls 2 notitle, \
  "../../data/is2/lndosst.dat"    u ($1/n):($2-$3) every 20 w p ls 3 notitle, \






unset multiplot
unset output
set terminal pop
reset
