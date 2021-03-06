set term pdfcairo size 40cm,30cm
set termopt enhanced
set lmargin 30
set rmargin 30
set bmargin 10
set tmargin 10
set key font ",30"
set xlabel font ",30"
set ylabel font ",30"
set tics font ",30"
set xtics offset 0,-1
set xlabel offset 0,-4
set ylabel offset -10,0
set key spacing 2
set ylabel "T*"
set xlabel "Timestep"
set output "baro_temp_compare.pdf"
plot "baro_temp_004_nochamal_1.dat" w l lt rgb "red" title "Shorter Run", "< tail -n 80000 baro_temp_004_new_1.dat" w l lt rgb "blue" title "Longer Run"
unset output

