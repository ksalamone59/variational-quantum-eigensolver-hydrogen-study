load "../style/default.gnu"
load "../style/term.gnu"
set out "time_scan.tex"
set lmargin 10
set rmargin 10

set title "Time Scan for VQE Optimization"
set xlabel '$r_{\text{max}}$'
set ylabel 'Number of Qubits' offset -2,0
set cblabel '$\log_{10}(\text{Average Time (s)})$' offset 4,0

set xrange[0.5:19.5]
set yrange[0.5:5.5]
# set cbrange[0:]

plot "time_scan.dat" using 2:1:(log10($3)) with image notitle