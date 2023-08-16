load '~/Projets/stage_M2/utils/setup.plot'

period1 = 100
period2 = 10
set xtics rotate by 45 right
set xlabel "step"

# ---------- Thermodynamic quantities ----------
set xlabel "step"

## Temperature
set ylabel "Température [K]"

set output 'relaxation_temp.pdf'
plot '../ch-sc/output/1/relaxation.log' every period2 using 1:2 with lines linestyle 1 title "chargé", \
     '../neutral/output/1/relaxation.log' every period2 using 1:2 with lines linestyle 2 title "neutre", \
     '../defected/output/1/relaxation.log' every period2 using 1:2 with lines linestyle 3 title "défectueux"

set output 'main_temp.pdf'
plot '../ch-sc/output/1/main.log' every period1 using 1:2 with lines linestyle 1 title "chargé", \
     '../neutral/output/1/main.log' every period1 using 1:2 with lines linestyle 2 title "neutre", \
     '../defected/output/1/main.log' every period1 using 1:2 with lines linestyle 3 title "défectueux"

## Potential energy
set ylabel "Énergie potentielle [kcal.mol^{-1}]"

set output 'relaxation_epot.pdf'
plot '../ch-sc/output/1/relaxation.log' every period2 using 1:3 with lines linestyle 1 title "chargé", \
     '../neutral/output/1/relaxation.log' every period2 using 1:3 with lines linestyle 2 title "neutre", \
     '../defected/output/1/relaxation.log' every period2 using 1:3 with lines linestyle 3 title "défectueux"

set output 'main_epot.pdf'
plot '../ch-sc/output/1/main.log' every period1 using 1:3 with lines linestyle 1 title "chargé", \
     '../neutral/output/1/main.log' every period1 using 1:3 with lines linestyle 2 title "neutre", \
     '../defected/output/1/main.log' every period1 using 1:3 with lines linestyle 3 title "défectueux"

## Pressure
set ylabel "Pression [atm]"

set output 'relaxation_press.pdf'
plot '../ch-sc/output/1/relaxation.log' every period2 using 1:5 with lines linestyle 1 title "chargé", \
     '../neutral/output/1/relaxation.log' every period2 using 1:5 with lines linestyle 2 title "neutre", \
     '../defected/output/1/relaxation.log' every period2 using 1:5 with lines linestyle 3 title "défectueux"

set output 'main_press.pdf'
plot '../ch-sc/output/1/main.log' every period1 using 1:5 with lines linestyle 1 title "chargé", \
     '../neutral/output/1/main.log' every period1 using 1:5 with lines linestyle 2 title "neutre", \
     '../defected/output/1/main.log' every period1 using 1:5 with lines linestyle 3 title "défectueux"

# ---------- Charges comparison ----------
set ylabel "<q> [e]"

## The charged and neutral
set output 'q_ch-sc-neutral.pdf'
plot '../ch-sc/output/1/q_lower.dat' every period1 using 1:2 with lines linestyle 1 title "chargé-", \
     '../ch-sc/output/1/q_upper.dat' every period1 using 1:2 with lines linestyle 2 title "chargé+", \
     '../neutral/output/1/q_upper.dat' every period1 using 1:2 with lines linestyle 4 title "neutre-", \
     '../neutral/output/1/q_lower.dat' every period1 using 1:2 with lines linestyle 5 title "neutre+"

set output 'q_ch-sc-neutral-out.pdf'
plot '../ch-sc/output/1/q_outer-lower.dat' every period1 using 1:2 with lines linestyle 1 title "chargé-", \
     '../ch-sc/output/1/q_outer-upper.dat' every period1 using 1:2 with lines linestyle 2 title "chargé+", \
     '../neutral/output/1/q_outer-upper.dat' every period1 using 1:2 with lines linestyle 4 title "neutre-", \
     '../neutral/output/1/q_outer-lower.dat' every period1 using 1:2 with lines linestyle 5 title "neutre+"

set output 'q_ch-sc-neutral-inn.pdf'
plot '../ch-sc/output/1/q_inner-lower.dat' every period1 using 1:2 with lines linestyle 1 title "chargé-", \
     '../ch-sc/output/1/q_inner-upper.dat' every period1 using 1:2 with lines linestyle 2 title "chargé+", \
     '../neutral/output/1/q_inner-upper.dat' every period1 using 1:2 with lines linestyle 4 title "neutre-", \
     '../neutral/output/1/q_inner-lower.dat' every period1 using 1:2 with lines linestyle 5 title "neutre+"

# ---------- RDFs comparison ----------
unset xtics
set xtics
set xlabel "r [Å]"
set ylabel "RDF"

## The charged and neutral
set output 'rdf_ch-sc-neutral-Na.pdf'
plot '../ch-sc/output/1/rdf_lower-Na.dat' using 1:2 with lines linestyle 1 title "chargé/négative", \
     '../ch-sc/output/1/rdf_upper-Na.dat' using 1:2 with lines linestyle 2 title "chargé/upper", \
     '../neutral/output/1/rdf_upper-Na.dat' using 1:2 with lines linestyle 4 title "neutre/négative", \
     '../neutral/output/1/rdf_lower-Na.dat' using 1:2 with lines linestyle 5 title "neutre/upper"

set output 'rdf_ch-sc-neutral-OH.pdf'
plot '../ch-sc/output/1/rdf_lower-OH.dat' using 1:2 with lines linestyle 1 title "chargé/négative", \
     '../ch-sc/output/1/rdf_upper-OH.dat' using 1:2 with lines linestyle 2 title "chargé/upper", \
     '../neutral/output/1/rdf_upper-OH.dat' using 1:2 with lines linestyle 4 title "neutre/négative", \
     '../neutral/output/1/rdf_lower-OH.dat' using 1:2 with lines linestyle 5 title "neutre/upper"

set output 'rdf_ch-sc-neutral-NaOH.pdf'
plot '../ch-sc/output/1/rdf_lower-Na.dat' using 1:2 with lines linestyle 1 title "chargé/Na-négative", \
     '../ch-sc/output/1/rdf_upper-OH.dat' using 1:2 with lines linestyle 2 title "chargé/OH-upper", \
     '../neutral/output/1/rdf_upper-Na.dat' using 1:2 with lines linestyle 4 title "neutre/Na-négative", \
     '../neutral/output/1/rdf_lower-OH.dat' using 1:2 with lines linestyle 5 title "neutre/OH-upper"
