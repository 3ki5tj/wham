#!/usr/bin/env gnuplot



set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced size 10, 7 font "Times, 24"
set output "whamcmp.eps"

set multiplot

# height of the bottom panel
ht0 = 0.2
# height of the middle panels
ht1 = 0.45
# height of the top panels
ht2 = 1 - ht0 - ht1

wd1 = 0.53
wd2 = 1 - wd1


set label 1 "(a)" at screen 0.01, 1 - 0.02                font "Times, 32"
set label 2 "(b)" at screen 0.01, ht1 + ht0 - 0.02        font "Times, 32"
set label 3 "(c)" at screen wd1 - 0.02, 1 - 0.02          font "Times, 32"
set label 4 "(d)" at screen wd1 - 0.02, ht1 + ht0 - 0.02  font "Times, 32"
set label 5 "(e)" at screen 0.01, ht0 - 0.0               font "Times, 32"


BOLTZ = 1.380658*6.0221367/1000

set size wd1, ht2
set origin 0.0, ht0 + ht1

set lmargin 7
set bmargin 0

set xtics 50 font "Times,20" offset 0, 0.5
set mxtics 5
set format x ""


#set logscale y
#set format y "10^{%T}"
set ytics 0.1 font "Times,20" offset 0.5, 0
set mytics 10
set ylabel "&{i}_{/Symbol*2.0 \341} {/Times-Italic f_i} - {/Times-Italic f_i}^{/*0.7 (ref)}&{i}_{/Symbol*2.0 \361}" offset 2.0, 0.0

color1  = "#224488"
color2  = "#aa2222"
color3  = "#804000"
color4  = "#aa00aa"
color5  = "#208020"
color6  = "#6688ff"
color7  = "#228888"

color9  = "#808080"

set style line 1  lw 2.0 lt 4 lc rgb color1  pt 5   ps 1.6
set style line 2  lw 2.0 lt 4 lc rgb color1  pt 13  ps 2.0
set style line 3  lw 2.0 lt 4 lc rgb color1  pt 7   ps 1.8
set style line 4  lw 4.0 lt 4 lc rgb color6  pt 3   ps 1.8

set style line 11 lw 2.0 lt 4 lc rgb color2  pt 9   ps 2.0
set style line 12 lw 2.0 lt 4 lc rgb color2  pt 11  ps 2.0
set style line 13 lw 2.0 lt 4 lc rgb color4  pt 15  ps 1.8
set style line 14 lw 2.0 lt 4 lc rgb color2  pt 13  ps 2.0

set style line 21 lw 3.0 lt 4 lc rgb color3  pt 1   ps 1.6
set style line 22 lw 3.0 lt 4 lc rgb color5  pt 2   ps 1.2
set style line 23 lw 3.0 lt 4 lc rgb color7  pt 3   ps 1.2
set style line 24 lw 3.0 lt 4 lc rgb color7  pt 15  ps 1.6

set style line 9  lw 2.0 lt 1 lc rgb color9  pt 10  ps 1.5

set xrange [305:415]

plot [:][-0.15:0.3] \
  "../../data/1VII/whamcmpr0.01/xvgerrwham.dat"          u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 1   notitle, \
  "../../data/1VII/whamcmpr0.01/xvgerrde100wham.dat"     u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 2   notitle, \
  "../../data/1VII/whamcmpr0.01/xvgerrmbar.dat"          u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 3   notitle, \
  "../../data/1VII/whamcmpr0.01/xvgerrbarmbar.dat"       u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 4   notitle, \
  "../../data/1VII/whamcmpr0.01/xvgerrstwham.dat"        u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 11  notitle, \
  "../../data/1VII/whamcmpr0.01/xvgerrde100stwham.dat"   u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 12  notitle, \
  "../../data/1VII/whamcmpr0.01/xvgerruiwham.dat"        u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 13  notitle, \
  "../../data/1VII/whamcmpr0.01/xvgerraveambar.dat"      u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 21  notitle, \
  "../../data/1VII/whamcmpr0.01/xvgerravebmbar.dat"      u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 22  notitle, \



set size wd2, ht2
set origin wd1, ht0 + ht1


set lmargin 3
unset ylabel


plot [:][-0.15:0.3] \
  "../../data/1VII/whamcmpr0.0001/xvgerrwham.dat"          u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 1   notitle, \
  "../../data/1VII/whamcmpr0.0001/xvgerrde100wham.dat"     u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 2   notitle, \
  "../../data/1VII/whamcmpr0.0001/xvgerrmbar.dat"          u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 3   notitle, \
  "../../data/1VII/whamcmpr0.0001/xvgerrbarmbar.dat"       u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 4   notitle, \
  "../../data/1VII/whamcmpr0.0001/xvgerrstwham.dat"        u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 11  notitle, \
  "../../data/1VII/whamcmpr0.0001/xvgerrde100stwham.dat"   u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 12  notitle, \
  "../../data/1VII/whamcmpr0.0001/xvgerruiwham.dat"        u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 13  notitle, \
  "../../data/1VII/whamcmpr0.0001/xvgerraveambar.dat"      u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 21  notitle, \
  "../../data/1VII/whamcmpr0.0001/xvgerravebmbar.dat"      u (1/BOLTZ/$2):(($1 > 0) ? ($3 -$5) : 1/0) w lp ls 22  notitle, \








#  Middle row





set size wd1, ht1
set origin 0, ht0

set tmargin 0.5
set bmargin 6
set lmargin 7

set xlabel "{/Times-Italic T_i} (K)" offset 1, 1
set format x "%g"

set ylabel "&{i}_{/Symbol*2.0 \341} {/Symbol-Oblique d} {/Times-Italic f_i}^2&{i}_{/Symbol*2.0 \361}^{1/2}" offset 3.0, 0.0
set ytics 0.05
set mytics 5

set key at 300, -0.020 left Left reverse font "Times,24" width -1 vertical maxrows 2 samplen 2.5 spacing 1.0


plot [:][0:0.1] \
  "../../data/1VII/whamcmpr0.01/xvgerrwham.dat"          u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 1   t "WHAM, {/Times-Italic h} = 1", \
  "../../data/1VII/whamcmpr0.01/xvgerrde100wham.dat"     u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 2   t "WHAM, {/Times-Italic h} = 100", \
  "../../data/1VII/whamcmpr0.01/xvgerrmbar.dat"          u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 3   t "MBAR", \
  "../../data/1VII/whamcmpr0.01/xvgerrbarmbar.dat"       u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 4   t "BAR, ({/Times-Italic i}, {/Times-Italic i} + 1)", \
  "../../data/1VII/whamcmpr0.01/xvgerrstwham.dat"        u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 11  t "ST-WHAM, {/Times-Italic h} = 1", \
  "../../data/1VII/whamcmpr0.01/xvgerrde100stwham.dat"   u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 12  t "ST-WHAM, {/Times-Italic h} = 100", \
  "../../data/1VII/whamcmpr0.01/xvgerruiwham.dat"        u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 13  t "UIM, {/Times-Italic h} = 1", \
  "../../data/1VII/whamcmpr0.01/xvgerraveambar.dat"      u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 21  t "Eq. (10)", \
  "../../data/1VII/whamcmpr0.01/xvgerravebmbar.dat"      u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 22  t "Eq. (11)", \


set size wd2, ht1
set origin wd1, ht0


set lmargin 3

unset ylabel
set ytics 0.5
set mytics 5


plot [:][0:1.0] \
  "../../data/1VII/whamcmpr0.0001/xvgerrwham.dat"          u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 1   notitle, \
  "../../data/1VII/whamcmpr0.0001/xvgerrde100wham.dat"     u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 2   notitle, \
  "../../data/1VII/whamcmpr0.0001/xvgerrmbar.dat"          u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 3   notitle, \
  "../../data/1VII/whamcmpr0.0001/xvgerrbarmbar.dat"       u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 4   notitle, \
  "../../data/1VII/whamcmpr0.0001/xvgerrstwham.dat"        u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 11  notitle, \
  "../../data/1VII/whamcmpr0.0001/xvgerrde100stwham.dat"   u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 12  notitle, \
  "../../data/1VII/whamcmpr0.0001/xvgerruiwham.dat"        u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 13  notitle, \
  "../../data/1VII/whamcmpr0.0001/xvgerraveambar.dat"      u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 21  notitle, \
  "../../data/1VII/whamcmpr0.0001/xvgerravebmbar.dat"      u (1/BOLTZ/$2):(($1 > 0) ? ($4) : 1/0) w lp ls 22  notitle, \






# Bottom row




reset

set size 1.0, ht0
set origin 0, 0

set tmargin 0
set lmargin 7

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
