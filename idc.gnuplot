set terminal png size 1400,1200
set output 'idc.png'
set nokey

#set logscale y
#set format y "10^{%L}"

#set yrange [1e-20:1]

plot 'idc.out' using 1:(abs($2)) ps 1.5 pt 6 w lp
