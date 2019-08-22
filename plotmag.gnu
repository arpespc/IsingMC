#!/usr/bin/gnuplot
#unset key
set term svg enhanced background rgb 'white'
set output "mag.svg"
plot './25by25/magnetization.dat' u 1:2 w lp pt 6 lw 2 title '25*25', \
     './20by20/magnetization.dat' u 1:2 w lp pt 6 lw 2 title '20*20', \
     './15by15/magnetization.dat' u 1:2 w lp pt 6 lw 2 title '15*15', \
     './12by12/magnetization.dat' u 1:2 w lp pt 6 lw 2 title '12*12', \
     './10by10/magnetization.dat' u 1:2 w lp pt 6 lw 2 title '10*10'
