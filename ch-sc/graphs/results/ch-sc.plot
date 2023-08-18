load '~/Projets/stage_M2/utils/setup.plot'

period = 100
period2 = 50
set xtics rotate by 45 right

# ---------- DZ ----------
set xlabel "step"
set ylabel "<dz> [Å]"

set output '1/ch-sc_dz-1.pdf'
plot '../../output/1/dz_lower.dat' every period using 1:2 with lines linestyle 1 title "négative", \
     '../../output/1/dz_upper.dat' every period using 1:2 with lines linestyle 2 title "upper"

set output '2/ch-sc_dz-1.pdf'
plot '../../output/2/dz_lower.dat' every period using 1:2 with lines linestyle 1 title "négative", \
     '../../output/2/dz_upper.dat' every period using 1:2 with lines linestyle 2 title "upper"

# ---------- Charges ----------
set ylabel "<q> [e]"

set output '1/ch-sc_q-1.pdf'
plot '../../output/1/q_outer-lower.dat' every period2 using 1:2 with lines linestyle 1 title "ext-", \
     '../../output/1/q_inner-lower.dat' every period2 using 1:2 with lines linestyle 4 title "int-", \
     '../../output/1/q_outer-upper.dat' every period2 using 1:2 with lines linestyle 2 title "ext+", \
     '../../output/1/q_inner-upper.dat' every period2 using 1:2 with lines linestyle 5 title "int+"

set output '2/ch-sc_q-2.pdf'
plot '../../output/2/q_outer-lower.dat' every period2 using 1:2 with lines linestyle 1 title "ext-", \
     '../../output/2/q_inner-lower.dat' every period2 using 1:2 with lines linestyle 4 title "int-", \
     '../../output/2/q_outer-upper.dat' every period2 using 1:2 with lines linestyle 2 title "ext+", \
     '../../output/2/q_inner-upper.dat' every period2 using 1:2 with lines linestyle 5 title "int+"

# ---------- Sodium density ----------
unset xtics
set xtics
set xlabel "z [Å]"
set ylabel "densité [Å^{-3}]"

stats '../../output/1/density_sodium.hist' using 1:2

set output 'ch-sc_density.pdf'
plot '../../output/1/density_sodium.hist' index 0 using ($2-STATS_min_y):3 with lines linestyle 1 title "départ", \
     '../../output/1/density_sodium.hist' index STATS_blocks-2 using ($2-STATS_min_y):3 with lines linestyle 2 title "fin"

stats '../../output/1/average_density_sodium.dat' using 1:2

set output 'ch-sc_average_density.pdf'
plot '../../output/1/average_density_sodium.dat' using ($1-STATS_min_x):2 with lines linestyle 1 notitle
