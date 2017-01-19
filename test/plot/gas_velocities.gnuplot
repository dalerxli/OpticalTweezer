folder="~/codes/OpticalTweezer/output/runs/170114_1401/"
set term pdfcairo size 40cm,30cm
set termopt enhanced
set xtics font ",30"
set ytics font ",30"
set xlabel font ",30"
set ylabel font ",30"
set key font ",30"
set xlabel "Timestep/100"
set ylabel "T*"
set key spacing 5
set bmargin 10
set lmargin 30
set tmargin 10
set rmargin 30
set xtics offset 0,-1
set xlabel offset 0,-3
set ylabel offset -6,0
set output "gas_temp.pdf"
plot folder."gasTempData.dat" using 1 w l lt rgb "red" lw 5 title "Incoming Particles", "" using 2 w l lt rgb "blue" lw 5 title "Outgoing Particles"
unset output
