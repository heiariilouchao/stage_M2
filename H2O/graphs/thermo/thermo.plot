load '~/Projets/stage_M2/utils/setup.plot'

main_period = 100
relaxation_period = 5
set xtics rotate by 45 right

# ---------- Temperature ----------
set xlabel "step"
set ylabel "Température [K]"

set output 'h2o_relaxation_temp.pdf'
plot 'reax_relaxation.log' every relaxation_period using 1:2 with lines linestyle 1 title "ReaxFF", \
     'spce_relaxation.log' every relaxation_period using 1:2 with lines linestyle 2 title "SPC/E"

set output 'h2o_main_temp.pdf'
plot 'reax.log' every main_period using 1:2 with lines linestyle 1 title "ReaxFF", \
     'spce.log' every main_period using 1:2 with lines linestyle 2 title "SPC/E"

# ---------- Pressure ----------
set ylabel "Pression [atm]"

set output 'h2o_relaxation_press.pdf'
plot 'reax_relaxation.log' every relaxation_period using 1:5 with lines linestyle 1 title "ReaxFF", \
     'spce_relaxation.log' every relaxation_period using 1:5 with lines linestyle 2 title "SPC/E"

set output 'h2o_main_press.pdf'
plot 'reax.log' every main_period using 1:5 with lines linestyle 1 title "ReaxFF", \
     'spce.log' every main_period using 1:5 with lines linestyle 2 title "SPC/E"

# ---------- Density ----------
set ylabel "Densité [g.cm^{-3}]"

set output 'h2o_relaxation_density.pdf'
plot 'reax_relaxation.log' every relaxation_period using 1:7 with lines linestyle 1 title "ReaxFF", \
     'spce_relaxation.log' every relaxation_period using 1:7 with lines linestyle 2 title "SPC/E"

# ---------- Potential Energy ----------
set ylabel "Énergie Potentielle (par molécule) [kcal.mol^{-1}]"

N_h2o = 267
stats 'reax_single.log' using 1:3

set output 'h2o_main_epot.pdf'
plot 'reax.log' every main_period using 1:($3 / N_h2o - STATS_mean_y) with lines linestyle 1 title "ReaxFF", \
     'spce.log' every main_period using 1:($3 / N_h2o) with lines linestyle 2 title "SPC/E"
