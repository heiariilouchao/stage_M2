set term pdfcairo font "Courrier, 11"
set xlabel "Step"


# The temperature
set ylabel "T [K]"
set output 'TEMP.pdf'
plot '../data/test.log' skip 1 u 1:2 with lines notitle


# The potential energy
set ylabel "E_p [kcal/mol]"
set output 'EPOT.pdf'
plot '../data/test.log' skip 1 u 1:3 with lines notitle

