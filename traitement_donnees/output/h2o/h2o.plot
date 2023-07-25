set term pdfcairo font "Cairo, 11pt"

set xlabel "décalage de configurations"
set ylabel "MSD [Å^2]"

set output 'MSD_raw.pdf'
plot 'reax-msd.dat' using 1:2 with lines title "ReaxFF", \
     'spce-msd.dat' using 1:2 with lines title "SPC/E"

set output 'MSD.pdf'
set xlabel "temps [s]"
set ylabel "MSD [m^2]"
plot 'reax-msd.dat' using ($1*1000*0.1e-15):($2*1e-20) with lines title "ReaxFF", \
     'spce-msd.dat' using ($1*1000*0.1e-15):($2*1e-20) with lines title "SPC/E"

set output 'MSD_annotated.pdf'
stats 'reax-msd.dat' using ($1*1000*0.1e-15):($2*1e-20) name 'REAX'
stats 'spce-msd.dat' using ($1*1000*0.1e-15):($2*1e-20) name 'SPCE'
plot 'reax-msd.dat' using ($1*1000*0.1e-15):($2*1e-20) with lines title sprintf("ReaxFF: %1.2e m^2.s^{-1}", REAX_slope), \
     'spce-msd.dat' using ($1*1000*0.1e-15):($2*1e-20) with lines title sprintf("SPC/E: %1.2e m^2.s^{-1}", SPCE_slope)

clear
set term pdfcairo font "Cairo, 11pt"

set xlabel "r [Å]"
set xrange [0:8]

set output 'h2o_rdf_oh.pdf'
set ylabel "g_{OH}"
plot 'reax-rdf.dat' using 1:3 with lines title "ReaxFF", \
     'spce-rdf.dat' using 1:3 with lines title "SPC/E", \
     'exp-rdf_oh.dat' using ($1+0.95):2 with lines title "Exp."

set output 'h2o_rdf_oo.pdf'
set ylabel "g_{OO}"
plot 'reax-rdf.dat' using 1:4 with lines title "ReaxFF", \
     'spce-rdf.dat' using 1:2 with lines title "SPC/E", \
     'exp-rdf_oo.dat' using ($1+0.88):2 with lines title "Exp."