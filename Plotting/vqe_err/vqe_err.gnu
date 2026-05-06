load "../style/default.gnu"
load "../style/term.gnu"
set out "vqe_err.tex"
set lmargin 10
set rmargin 10

set title "Variational Error"
set xlabel '$r_{\text{max}}$'
set ylabel 'Number of Qubits' offset -2,0
set cblabel '$\log_{10}(\text{Error})$' offset 4,0

set xrange[0.5:19.5]
set yrange[0.5:5.5]

# set logscale cb
set format cb "%.1e"

plot "vqe_err.dat" using 2:1:3 with image notitle