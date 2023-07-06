set term pdfcairo font "Courrier, 11"


# Radial distribution functions
set output "RDF.pdf"
set xlabel "r [Å]"
set ylabel "RDF"
plot '../data/rdf.dat' using 1:2 with lines title "g_{OH}", \
     '../data/rdf.dat' using 1:3 with lines title "g_{OO}", \
     '../data/rdf.dat' using 1:4 with lines title "g_{OH}"

set output "RDF-comp.pdf"
plot '../data/rdf.dat' using 1:2 with lines title "g_{OH}", \
     '../data/rdf.dat' using 1:3 with lines title "g_{OO}", \
     '../data/rdf.dat' using 1:4 with lines title "g_{OH}", \
     '../data/RDF-OH.dat' using 1:2 with lines title "g_{OH} (vmd)", \
     '../data/RDF-OO.dat' using 1:2 with lines title "g_{OO} (vmd)", \
     '../data/RDF-HH.dat' using 1:2 with lines title "g_{HH} (vmd)"

set output "RDF-comp-shifted.pdf"
plot '../data/rdf.dat' using ($1-0.1):2 with lines title "g_{OH}", \
     '../data/rdf.dat' using ($1-0.1):3 with lines title "g_{OO}", \
     '../data/rdf.dat' using ($1-0.1):4 with lines title "g_{OH}", \
     '../data/RDF-OH.dat' using 1:2 with lines title "g_{OH} (vmd)", \
     '../data/RDF-OO.dat' using 1:2 with lines title "g_{OO} (vmd)", \
     '../data/RDF-HH.dat' using 1:2 with lines title "g_{HH} (vmd)"


# Mean Squared Displacements
set output "MSD.pdf"
set xlabel "step"
set ylabel "MSD [Å²]"
plot '../data/msd.dat' using 1:2 with lines notitle

set output "RMSD.pdf"
set ylabel "RMSD [Å]"
plot '../data/msd.dat' using 1:3 with lines notitle

set output "RMSD-comp.pdf"
plot '../data/msd.dat' using 1:3 with lines title "custom", \
     '../data/RMSD.dat' using 1:2 with lines title "vmd"

