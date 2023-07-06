# T = 300 K
# P = 0 atm
# t_{eq} =! 5e-10 s / 5000000 steps


# Graphical properties
set key font "Cairo,14"
set linetype 1 lw 2
set linetype 2 lw 2


# MSD
set title "Mean Squared Displacement" font "Cairo,14"
set xlabel "{/Symbol t} [s]" font "Cairo,14"
set ylabel "MSD [m^2]" font "Cairo,14"
stats [5.e-13:] '../data/spce/spce.msd' using ($1*1000.0*0.1e-15):($2*1.e-20)
stats [5.e-13:] '../data/reaxff/reax.msd' using ($1*1000.0*0.1e-15):($2*1.e-20)
plot '../data/spce/spce.msd' using ($1*1000.0*0.1e-15):($2*1.e-20) with lines title "SPC/E", \
     '../data/reaxff/reax.msd' using ($1*1000.0*0.1e-15):($2*1.e-20) with lines title "ReaxFF"
pause -1


# RDF
set title "Radial Distribution Function: O-O" font "Cairo,14"
set xlabel "r [Å]" font "Cairo,14"
set ylabel "g" font "Cairo,14"
plot '../data/spce/spce.rdf' using 1:4 with lines title "SPC/E", \
     '../data/reaxff/reax.rdf' using 1:4 with lines title "ReaxFF"
pause -1


set title "Radial Distribution Function: O-H" font "Cairo,14"
plot '../data/spce/spce.rdf' using 1:3 with lines title "SPC/E", \
     '../data/reaxff/reax.rdf' using 1:3 with lines title "ReaxFF"
pause -1

