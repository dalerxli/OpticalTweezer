set term pdfcairo size 40cm,30cm
set yrange [*:*]
set xrange [*:*]
set termopt enhanced
set xlabel "Time"
set ylabel "T*"
set y2label "E_{tot}"
set tics font ",30"
set xlabel font ",30"
set ylabel font ",30"
set y2label font ",30"
set bmargin 15
set lmargin 20
set tmargin 10
set rmargin 20 
set xlabel offset 0,-7
set ylabel offset -7,0
set y2label offset 7,0
set y2tics
set ytics nomirror
set xtics offset 0,-3
set key left
set key font ",30"
set key spacing 5
set key at 1000,0.89
set output "ehex_te.pdf"
#plot "ehex_temp_004_1.dat" using ($0*0.0217):1 w l axes x1y1 lw 2 lt rgb "red" title "T*" ,\
#     "ehex_energy_004_1.dat" using ($0*0.0217):1 w l axes x1y2 lw 5 lt rgb "blue" title "E_{tot}" 
plot "ehex_temp_000_new_1_new.dat" using ($0*0.01):1 w l axes x1y1 lw 2 lt rgb "red" title "T*" ,\
     "ehex_energy_000_new_1_new.dat" using ($0*0.01):1 w l axes x1y2 lw 5 lt rgb "blue" title "E_{tot}" 
unset output
