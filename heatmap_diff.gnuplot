set terminal png size 1400,1200
set output 'heatmap_diff.png'
set nokey

#set palette negative

set logscale cb
set format cb "10^{%L}"

#load 'blackbody.pal'
set palette cubehelix

#### bump parameters.
#set xrange [-4:4]
#set yrange [-4:4]

#set cbrange [1e0:1e-4]

#### lattice parameters.
set xrange [0:7]
set yrange [0:7]

set cbrange [1e-7:1e2]

#### hohlraum parameters.
# set xrange [0:1.3]
# set yrange [1.3:0]

# set cbrange [1e0:1e-6]
#set cbrange [0:100]

set arrow from 0,0.0625 to 7,0.0625 lc "white" front nohead
set arrow from 0,2.375  to 7,2.375  lc "white" front nohead
set arrow from 0,4.6875 to 7,4.6875 lc "white" front nohead

plot 'heatmap.out' using 1:2:(abs($3)) with image
