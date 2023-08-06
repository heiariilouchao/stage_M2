reset

set term pdfcairo font "Cairo,14"

set linetype 1 linecolor rgb "red" linewidth 2
set linetype 2 linecolor rgb "blue" linewidth 2
set linetype 3 linecolor rgb "orange" linewidth 2

set style line 1 linetype 1
set style line 2 linetype 2
set style line 3 linetype 3
set style line 4 linetype 1 dashtype 2
set style line 5 linetype 2 dashtype 2
set style line 6 linetype 3 dashtype '..'

set grid
set key outside
