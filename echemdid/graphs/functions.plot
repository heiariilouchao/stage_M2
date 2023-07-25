set term pdfcairo font "Cairo,14"

set linetype 1 linecolor rgb "red" linewidth 2
set linetype 2 linecolor rgb "blue" linewidth 2
set grid

# ---------- w ----------
Rc = 4.0
print sprintf("Rc = %f", Rc)
# max_W = 3.0 * w(1.42) + 6 * w(2.46) + 3 * w(2.83) + 6 * w(3.75) + 2 / 3 * (w(3.42) + 2 * w(3.57) + 2 * w(3.84) + 2 * w(3.97))
max_W = 3.0
print sprintf("max_W = %f", max_W)
N = 2.0 * 3 * 540 / (540 * max_W)
print sprintf("N = %f", N)
w(x) = (x < Rc) ? N * (1.0 - (x / Rc)**2)**2 : 0.0

set xlabel "R"
set ylabel "w(R)"
set output 'echemdid_w.pdf'
plot [0:Rc] w(x) with lines linetype 1 notitle


# ---------- F ----------
W0 = w(0.99 * Rc)
print sprintf("W0 = %f", W0)
F(x) = (x < W0) ? (1.0 - (x / W0)**2)**2 : 0.0

set xlabel "W"
set ylabel "F(W)"
set output 'echemdid_F.pdf'
plot [max_W:0] F(x) with lines linetype 2 notitle


# ---------- w and F ----------
set xlabel "R"
unset ylabel
set output 'echemdid_w_F.pdf'
plot [0:Rc] w(x) with lines linetype 1 title "w(R)", \
            F(w(x)) with lines linetype 2 title "F(R)"
