load '~/Projets/stage_M2/utils/setup.plot'
period = 100
period2 = 20
set xtics rotate by 45 right

## Temperature
set xlabel "t"
set ylabel "Température [K]"

set output 'defected_relaxation_temp.pdf'
plot '../../output/1/relaxation.log' every period2 using 1:2 with lines linestyle 1 title "1", \
     '../../output/2/relaxation.log' every period2 using 1:2 with lines linestyle 2 title "2"

set output 'defected_main_temp.pdf'
plot '../../output/1/main.log' every period using 1:2 with lines linestyle 1 title "1", \
     '../../output/2/main.log' every period using 1:2 with lines linestyle 2 title "2"

## Pressure
set ylabel "Pression [atm]"

set output 'defected_relaxation_press.pdf'
plot '../../output/1/relaxation.log' every period2 using 1:5 with lines linestyle 1 title "1", \
     '../../output/2/relaxation.log' every period2 using 1:5 with lines linestyle 2 title "2"

set output 'defected_main_press.pdf'
plot '../../output/1/main.log' every period using 1:5 with lines linestyle 1 title "1", \
     '../../output/2/main.log' every period using 1:5 with lines linestyle 2 title "2"

## Potential Energy
set ylabel "Énergie potentielle totale[kcal.mol^{-1}]"

set output 'defected_relaxation_epot.pdf'
plot '../../output/1/relaxation.log' every period2 using 1:3 with lines linestyle 1 title "1", \
     '../../output/2/relaxation.log' every period2 using 1:3 with lines linestyle 2 title "2"

set output 'defected_main_epot.pdf'
plot '../../output/1/main.log' every period using 1:3 with lines linestyle 1 title "1", \
     '../../output/2/main.log' every period using 1:3 with lines linestyle 2 title "2"

## Density
set ylabel "Densité [g.cm^{-3}]"

set output 'defected_relaxation_density.pdf'
plot '../../output/1/relaxation.log' every period2 using 1:7 with lines linestyle 1 title "1", \
     '../../output/2/relaxation.log' every period2 using 1:7 with lines linestyle 2 title "2"
