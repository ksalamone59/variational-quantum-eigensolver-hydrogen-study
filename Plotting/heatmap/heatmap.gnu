load "../style/default.gnu"
load "../style/term.gnu"
set out "heatmap.tex"

set lmargin 10
set rmargin 7

set xlabel '$r_{\text{max}}$'
set ylabel 'Number of Qubits' offset -2,0
set cblabel '$\log_{10}(\text{Error})$' offset 4,0
set title 'Error Heatmap for VQE Simulation'

set xrange[0.5:9.5]
set yrange[0.5:19.5]

# ytics from 1 to 19 by 1
set ytics 1, 2, 19

set key opaque
set key box

plot "heatmap.dat" using 1:2:3 with image notitle,\
    '-' using 1:2 with p ls 2 title "Minimum Error: $4.29e-5$    "
    9 9
    e