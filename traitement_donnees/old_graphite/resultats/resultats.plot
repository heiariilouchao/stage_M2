
# Averaged
set title "Average charges comparison"
set xlabel "iterations"
set ylabel "q/N [e]"

plot '../graphite.hist' using 1:2 with lines title "out+", \
     '' using 1:4 with lines title "out-", \
     '' using 1:6 with lines title "inn+", \
     '' using 1:8 with lines title "inn-", \
     '' using 1:10 with lines title "b = 2"
pause -1


# Samples
set title "Sample charges"
set xlabel "iterations"
set ylabel "q [e]"

plot '../samples.log' using 1:2 with lines title "out+", \
     '' using 1:3 with lines title "out-", \
     '' using 1:4 with lines title "inn+", \
     '' using 1:5 with lines title "inn-", \
     '' using 1:6 with lines title "b = 2"
pause -1

plot '../samples.log' using 1:2 with lines title "out+", \
     '' using 1:4 with lines title "inn+", \
     '' using 1:6 with lines title "b = 2"
pause -1

plot '../samples.log' using 1:3 with lines title "out-", \
     '' using 1:5 with lines title "inn-"
pause -1
