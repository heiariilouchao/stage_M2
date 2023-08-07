load '~/Projets/stage_M2/utils/setup.plot'

period = 100
period2 = 50

# ---------- DZ ----------
set xlabel "step"
set ylabel "<dz> [Å]"

set output 'ch-sc_dz.pdf'
plot 'dz_negative.dat' every period using 1:2 with lines linestyle 1 title "négative", \
     'dz_positive.dat' every period using 1:2 with lines linestyle 2 title "positive"

# ---------- Charges ----------
set ylabel "<q> [e]"

set output 'ch-sc_q.pdf'
plot 'q_outer-negative.dat' every period2 using 1:2 with lines linestyle 1 title "ext-", \
     'q_inner-negative.dat' every period2 using 1:2 with lines linestyle 2 title "int-", \
     'q_outer-positive.dat' every period2 using 1:2 with lines linestyle 4 title "ext+", \
     'q_inner-positive.dat' every period2 using 1:2 with lines linestyle 5 title "int+"
