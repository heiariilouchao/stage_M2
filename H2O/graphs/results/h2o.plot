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


# ---------- MSDs ----------
coeff_t = 1000.0*1.0e-15
coeff_MSD = 1.0e-20

set xtics rotate by 30 right
set xlabel "temps [s]"
set ylabel "MSD [m^2]"
set output 'h2o_msd.pdf'
stats 'reax-msd.dat' using ($1 * coeff_t):($2 * coeff_MSD) name 'REAX'
stats 'spce-msd.dat' using ($1 * coeff_t):($2 * coeff_MSD) name 'SPCE'
plot 'reax-msd.dat' using ($1 * coeff_t):($2 * coeff_MSD) with lines linestyle 1 title "ReaxFF", \
     'spce-msd.dat' using ($1 * coeff_t):($2 * coeff_MSD) with lines linestyle 2 title "SPC/E"


# ---------- RDFs ----------
rc = 8.0

set key inside
set xtics rotate by 0 center
set xlabel "r [Å]"
set xrange [0:8]

set output 'h2o_rdf_oh.pdf'
set ylabel "g_{OH}"
plot 'reax-rdf.dat' using 1:3 with lines title "ReaxFF", \
     'spce-rdf.dat' using 1:3 with lines title "SPC/E", \
     'exp-rdf_oh.csv' using 1:2 with lines title "Exp."

set output 'h2o_rdf_oo.pdf'
set ylabel "g_{OO}"
plot 'reax-rdf.dat' using 1:4 with lines linestyle 1 title "ReaxFF", \
     'spce-rdf.dat' using 1:2 with lines linestyle 2 title "SPC/E", \
     'exp-rdf_oo.csv' using 1:2 with lines linestyle 3 title "Exp."