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
set yrange [10:4000]

set logscale y2
set y2tics nomirror font "Times,20" offset -0.5, 0
set y2label "Relative run time, {/Times-Italic t}_{/=16 rel} = {/Times-Italic t}/{/Symbol-Oblique t}" textcolor rgb "#404040" offset -3.5, 0
set y2range [0.01:4]
set rmargin 5.0

set key invert font "Times,24" spacing 1.2
#set key Left reverse spacing 1.5

# Color figure

#color1  = "#224488"
#color2  = "#aa0000"
#color3  = "#aa00aa"
#color4  = "#00aa00"
#color5  = "#00aaaa"

#color11 = "#aaaaff"
#color12 = "#ffaaaa"
#color13 = "#ffaaff"
#color14 = "#aaffaa"
#color15 = "#aaffff"

#set style line 1  lw 2.0 lt 1 lc rgb color1  pt 4   ps 1.8
#set style line 2  lw 2.0 lt 1 lc rgb color2  pt 8   ps 2.0
#set style line 3  lw 2.0 lt 1 lc rgb color3  pt 10  ps 2.0
#set style line 4  lw 2.0 lt 1 lc rgb color4  pt 6   ps 2.0
#set style line 5  lw 2.0 lt 1 lc rgb color5  pt 12  ps 2.0
#
#set style line 11 lw 2.0 lt 1 lc rgb color11 pt 5   ps 1.8
#set style line 12 lw 2.0 lt 1 lc rgb color12 pt 9   ps 2.0
#set style line 13 lw 2.0 lt 1 lc rgb color13 pt 11  ps 2.0
#set style line 14 lw 2.0 lt 1 lc rgb color14 pt 7   ps 2.0
#set style line 15 lw 2.0 lt 1 lc rgb color15 pt 13  ps 2.0


# black and white figures

color1  = "#000000"
color2  = "#111111"
color3  = "#222222"
color4  = "#333333"
color5  = "#444444"

color11 = "#888888"
color12 = "#aaaaaa"
color13 = "#bbbbbb"
color14 = "#cccccc"
color15 = "#dddddd"

set style line 1  lw 2.0 lt 1 lc rgb color1  pt 4   ps 1.8
set style line 2  lw 2.0 lt 2 lc rgb color2  pt 8   ps 2.0
set style line 3  lw 2.0 lt 5 lc rgb color3  pt 10  ps 2.0
set style line 4  lw 2.0 lt 4 lc rgb color4  pt 6   ps 2.0
set style line 5  lw 2.0 lt 5 lc rgb color5  pt 12  ps 2.0

set style line 11 lw 2.0 lt 1 lc rgb color11 pt 5   ps 1.8
set style line 12 lw 2.0 lt 2 lc rgb color12 pt 9   ps 2.0
set style line 13 lw 2.0 lt 5 lc rgb color13 pt 11  ps 2.0
set style line 14 lw 2.0 lt 4 lc rgb color14 pt 7   ps 2.0
set style line 15 lw 2.0 lt 5 lc rgb color15 pt 13  ps 2.0

is2it = real(`sed '2!d' ../../data/is2/run/is2wham.dat | awk '{print $2}'`)
is2tm = real(`sed '2!d' ../../data/is2/run/is2tmwham.dat | awk '{print $2}'`)
is2scl = is2tm*1000.0/is2it
print "Ising scaling: ", is2scl

vlwit = real(`sed '2!d' ../../data/1VII/xvgwham.dat | awk '{print $2}'`)
vlwtm = real(`sed '2!d' ../../data/1VII/xvgtmwham.dat | awk '{print $2}'`)
vlwscl = vlwtm*1000.0/vlwit
print "Villin WHAM scaling: ", vlwscl

vlmit = real(`sed '2!d' ../../data/1VII/xvgmbar.dat | awk '{print $2}'`)
vlmtm = real(`sed '2!d' ../../data/1VII/xvgtmmbar.dat | awk '{print $2}'`)
vlmscl = vlmtm*1000.0/vlmit
print "Villin MBAR scaling: ", vlmscl

lj2it = real(`sed '2!d' ../../data/lj/lj2wham.dat | awk '{print $2}'`)
lj2tm = real(`sed '2!d' ../../data/lj/lj2tmwham.dat | awk '{print $2}'`)
lj2scl = lj2tm*1000.0/lj2it
print "LJ2 scaling: ", lj2scl

plot [0:20][] \
  "../../data/1VII/xvgtmmbar.dat"     u 1:(($1 > -1 && $1 <= 12) ? ($2/vlmscl) : 1/0) axes x1y2 w lp ls 13 t "MBAR, villin, {/Times-Italic NVT}, {/Times-Italic t}_{/=16 rel}", \
  "../../data/1VII/xvgtmwham.dat"     u 1:(($1 > -1 && $1 <= 12) ? ($2/vlwscl) : 1/0) axes x1y2 w lp ls 12 t "WHAM, villin, {/Times-Italic NVT}, {/Times-Italic t}_{/=16 rel}", \
  "../../data/lj/lj2tmwham.dat"       u 1:(($1 > -1) ? ($2/lj2scl) : 1/0) axes x1y2 w lp ls 14 t "WHAM, LJ, {/Times-Italic NpT}, {/Times-Italic t}_{/=16 rel}", \
  "../../data/is2/run/is2tmwham.dat"  u 1:(($1 > -1) ? ($2/is2scl) : 1/0) axes x1y2 w lp ls 11 t "WHAM, Ising, {/Times-Italic t}_{/=16 rel}", \
  "../../data/1VII/xvgmbar.dat"       u 1:(($1 > -1 && $1 <= 12) ? ($2) : 1/0) w lp ls 3 t "MBAR, villin, {/Times-Italic NVT}, {/Times-Italic N}@_{/=16 iter}   ", \
  "../../data/1VII/xvgwham.dat"       u 1:(($1 > -1 && $1 <= 12) ? ($2) : 1/0) w lp ls 2 t "WHAM, villin, {/Times-Italic NVT}, {/Times-Italic N}@_{/=16 iter}   ", \
  "../../data/lj/lj2wham.dat"         u 1:(($1 > -1) ? ($2) : 1/0) w lp ls 4 t "WHAM, LJ, {/Times-Italic NpT}, {/Times-Italic N}@_{/=16 iter}   ", \
  "../../data/is2/run/is2wham.dat"    u 1:(($1 > -1) ? ($2) : 1/0) w lp ls 1 t "WHAM, Ising, {/Times-Italic N}@_{/=16 iter}   ", \



unset output
set terminal pop
reset
