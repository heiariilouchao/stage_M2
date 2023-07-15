set term gif animate delay 4
set output 'sodium.gif'
stats '../sodium.hist' u 2:3
set xlabel "z [Å]"
set ylabel "density [Å^{-3}]"
set xrange [STATS_min_x:STATS_max_x]
set yrange [STATS_min_y:STATS_max_y]
do for [i=0:int(STATS_blocks-2):int(STATS_blocks/600)] \
{
    plot '../sodium.hist' index i using 2:3 with lines title sprintf("step: %f", $1)
}