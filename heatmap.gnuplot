set terminal png size 1400,1200
set output 'heatmap.png'
set nokey

#set palette negative

set logscale cb
set format cb "10^{%L}"

#load 'jet.pal'
#load 'blackbody.pal'
#load 'viridis.pal'

set palette cubehelix

#### bump parameters.
# set xrange [0:4]
# set yrange [0:4]

# set cbrange [1e1:1e-3]

#### lattice parameters.
set xrange [0:7]
set yrange [7:0]

set cbrange [1e-6:1e0]

#### hohlraum parameters.
#set xrange [0:1.3]
#set yrange [1.3:0]

#set cbrange [1e0:1e-4]
#set cbrange [0:100]

plot 'heatmap.out' using 1:2:(abs($3)) with image
