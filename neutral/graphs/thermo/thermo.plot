load '~/Projets/stage_M2/utils/setup.plot'
period = 100
period2 = 20
set xtics rotate by 45 right

## Temperature
set xlabel "t"
set ylabel "Température [K]"

set output 'neutral_relaxation_temp.pdf'
plot '../../output/relaxation1.log' every period2 using 1:2 with lines linestyle 1 title "1", \
     '../../output/relaxation2.log' every period2 using 1:2 with lines linestyle 2 title "2"

set output 'neutral_main_temp.pdf'
plot '../../output/main1.log' every period using 1:2 with lines linestyle 1 title "1", \
     '../../output/main2.log' every period using 1:2 with lines linestyle 2 title "2"

## Pressure
set ylabel "Pression [atm]"

set output 'neutral_relaxation_press.pdf'
plot '../../output/relaxation1.log' every period2 using 1:5 with lines linestyle 1 title "1", \
     '../../output/relaxation2.log' every period2 using 1:5 with lines linestyle 2 title "2"

set output 'neutral_main_press.pdf'
plot '../../output/main1.log' every period using 1:5 with lines linestyle 1 title "1", \
     '../../output/main2.log' every period using 1:5 with lines linestyle 2 title "2"

## Potential Energy
set ylabel "Énergie potentielle totale [kcal.mol^{-1}]"

set output 'neutral_relaxation_epot.pdf'
plot '../../output/relaxation1.log' every period2 using 1:3 with lines linestyle 1 title "1", \
     '../../output/relaxation2.log' every period2 using 1:3 with lines linestyle 2 title "2"

set output 'neutral_main_epot.pdf'
plot '../../output/main1.log' every period using 1:3 with lines linestyle 1 title "1", \
     '../../output/main2.log' every period using 1:3 with lines linestyle 2 title "2"

## Density
set ylabel "Densité [g.cm^{-3}]"

set output 'neutral_relaxation_density.pdf'
plot '../../output/relaxation1.log' every period2 using 1:7 with lines linestyle 1 title "1", \
     '../../output/relaxation2.log' every period2 using 1:7 with lines linestyle 2 title "2"
