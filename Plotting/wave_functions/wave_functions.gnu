load "../style/default.gnu"
load "../style/term.gnu"
set out "wave_functions.tex"
set title '$u(r)$ Wave Functions from Classical Approach'
set xlabel 'r (m)'
set ylabel '$u(r)$ (normalized)' 
a0 = 5.29177e-11
gs(x) = 2.0 * (a0**(-1.5)) * x * exp(-x/a0)
excited_state(x) = (1/(2*sqrt(2))) * (a0**(-1.5))*x*(2-x/a0)*exp(-x/(2*a0))
set xrange[0:1e-9]
set format y "%1.e"
plot gs(x) w p ls 2 title "Analytical Ground State",\
    excited_state(x) w p ls 4 title "Analytical First Excited State",\
    "../../data/wavefunctions.dat" u 1:2 w l ls 1 title "Ground State",\
    "../../data/wavefunctions.dat" u 1:3 w l ls 3 title "First Excited State"