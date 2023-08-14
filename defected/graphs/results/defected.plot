load '~/Projets/stage_M2/utils/setup.plot'

period = 100
period2 = 50
set xtics rotate by 45 right

# ---------- DZ ----------
set xlabel "step"
set ylabel "<dz> [Å]"

set output '1/defected_dz.pdf'
plot '../../output/1/dz_negative.dat' every period using 1:2 with lines linestyle 1 title "négative", \
     '../../output/1/dz_positive.dat' every period using 1:2 with lines linestyle 2 title "positive"

set output '2/defected_dz.pdf'
plot '../../output/2/dz_negative.dat' every period using 1:2 with lines linestyle 1 title "négative", \
     '../../output/2/dz_positive.dat' every period using 1:2 with lines linestyle 2 title "positive"

# ---------- Charges ----------
set ylabel "<q> [e]"

set output '1/defected_q.pdf'
plot '../../output/1/q_outer-negative.dat' every period2 using 1:2 with lines linestyle 1 title "ext-", \
     '../../output/1/q_inner-negative.dat' every period2 using 1:2 with lines linestyle 4 title "int-", \
     '../../output/1/q_outer-positive.dat' every period2 using 1:2 with lines linestyle 2 title "ext+", \
     '../../output/1/q_inner-positive.dat' every period2 using 1:2 with lines linestyle 5 title "int+"

set output '2/defected_q.pdf'
plot '../../output/2/q_outer-negative.dat' every period2 using 1:2 with lines linestyle 1 title "ext-", \
     '../../output/2/q_inner-negative.dat' every period2 using 1:2 with lines linestyle 4 title "int-", \
     '../../output/2/q_outer-positive.dat' every period2 using 1:2 with lines linestyle 2 title "ext+", \
     '../../output/2/q_inner-positive.dat' every period2 using 1:2 with lines linestyle 5 title "int+"

# ---------- Sodium density ----------
set xlabel "z [Å]"
set ylabel "densité [Å^{-3}]"

stats '../../output/1/density_sodium.hist' using 1:2

set output '1/defected_sodium-density.pdf'
plot '../../output/1/density_sodium.hist' index 0 using 2:3 with lines linestyle 1 title "départ", \
     '../../output/1/density_sodium.hist' index STATS_blocks-2 using 2:3 with lines linestyle 2 title "fin"

stats '../../output/2/density_sodium.hist' using 1:2

set output '2/defected_sodium-density.pdf'
plot '../../output/2/density_sodium.hist' index 0 using 2:3 with lines linestyle 1 title "départ", \
     '../../output/2/density_sodium.hist' index STATS_blocks-2 using 2:3 with lines linestyle 2 title "fin"
