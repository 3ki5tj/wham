#!/usr/bin/env gnuplot



set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced size 10, 4.0 font "Times, 40"
set output "gausconv.eps"
set multiplot


wl = 0.45
wr = 1 - wl


n = 64*64



ldx = 0.013
ldy = 0.06
set label 10 "(a)" at screen ldx,    1 - ldy font "Times, 48"
set label 11 "(b)" at screen wl+ldx, 1 - ldy font "Times, 48"


# for color figure

color1  = "#aa0000"
color2  = "#224488"
color3  = "#006000"
color4  = "#aa00aa"

# for black/white figure

#color1 = "#000000"
#color2 = "#222222"
#color3 = "#444444"
#color4 = "#666666"

color9  = "#202020"

set style line 1  lw 4.0 lt 1 lc rgb color1  pt 10  ps 1.6
set style line 2  lw 4.0 lt 2 lc rgb color2  pt 8   ps 1.6
set style line 3  lw 4.0 lt 3 lc rgb color3  pt 4   ps 1.2
set style line 4  lw 4.0 lt 5 lc rgb color4  pt 6   ps 1.2

set style line 9  lw 4.0 lt 1 lc rgb color9  pt 1   ps 1.0

set size wl, 0.9
set origin 0, 0


set lmargin 0.1
set rmargin 0.0
set tmargin 1.0
set bmargin 2.8


unset border

set xlabel "{/Times-Italic E}"
#unset xlabel
unset xtics

unset logscale
unset ylabel
unset ytics

#set key Left reverse


ec = 0.5;
xmax = 4;
ymax = 0.45;

set label 1 "{/Symbol-Oblique r}_1({/Times-Italic E})" at  ec+0.8, 0.35 front
set label 2 "{/Symbol-Oblique r}_2({/Times-Italic E})" at -ec-2.0, 0.35 front
set label 3 "{/Times-Italic A}_{12}" at -0.1, 0.07 front

set arrow 1 from 0,0 to 0,ymax*1.2 back ls 9 nohead
set arrow 2 from -xmax,0 to xmax+0.3,0 back ls 9 nohead

plot [-xmax:xmax][0:ymax] \
  1/sqrt(2*pi)*exp(-(x-ec)**2/2) ls 1 notitle, \
  1/sqrt(2*pi)*exp(-(x+ec)**2/2) ls 2 notitle, \
  1/sqrt(2*pi)*exp(-x*x/2-ec*ec/2)/2/cosh(x*ec) w filledcurves lc rgb "#e0e0e0" notitle


unset label 1
unset label 2
unset label 3
unset arrow 1
unset arrow 2



set size wr, 1
set origin wl, 0

set lmargin 7
set rmargin 1

set border

set xlabel "{/Symbol-Oblique s}_{/Symbol-Oblique b} {/Symbol-Oblique s}_{/Times-Italic E}" offset 0, 0.5
set xtics 1.0 offset 0, 0.3
set mxtics 2

set logscale y
set format y "10^{%T}"
set ylabel "1 - {/Symbol-Oblique l}_1" offset 0.5, 0.0
set ytics 10 offset 0.5, 0
set mytics 10

set key at -0.05, 2e-6 left bottom Left reverse


rr(b, t) = exp(-b*b/2)/(1+0.5*exp(t*b*b))/2/sqrt(1 + b*b/(1+0.5*exp(t*b*b)))*(1 + b**4*(2-0.5*exp(t*b*b))/8/(1+b*b+0.5*exp(t*b*b)))

plot [0:5][1e-6:] \
  exp(-x*x/2)/sqrt(1+x*x)*(1+x**4/4/(1+x*x)**2) ls 3 t "Two-temperature", \
  rr(sqrt(3.0/8)*x,-2.5) + 2*rr(sqrt(1.5)*x,0.5) ls 2 t "Three-temperature", \
  1/(1 + x * x) ls 1 t "Continuous-temperature"


unset multiplot
unset output
set terminal pop
reset
