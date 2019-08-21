#!/usr/bin/gnuplot

set term svg enhanced background rgb 'white'
set output "spin.svg"
# wxt
#set terminal wxt size 350,262 enhanced font 'Verdana,10' persist
# png
#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'
#set output 'heat_map_interpolation1.png'

set border linewidth 0
unset key
#unset colorbox
#unset tics
#set lmargin screen 0.1
#set rmargin screen 0.9
#set tmargin screen 0.9
#set bmargin screen 0.1
set cbrange [-6:6]
#set palette grey

set pm3d map
splot './M_mat_10.000000' matrix
