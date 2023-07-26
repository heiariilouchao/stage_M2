set term pdfcairo font "Cairo,14"

set linetype 1 linecolor rgb "red" linewidth 2
set linetype 4 linecolor rgb "black"
set grid


# ---------- SPC/E's Lennard-Jones ----------
A = 0.37122
B = 0.3428
Rc = 1.0
min_R = 0.3

U(x) = -(A / x)**6 + (B / x)**12

set xlabel "r [nm]"
set ylabel "U [kcal.mol^{-1}]"
set output 'spce_lj.pdf'
plot [min_R:Rc] U(x) with lines linetype 1 notitle, \
                0.0 with lines linetype 4 notitle