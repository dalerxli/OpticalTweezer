set term pdfcairo size 40cm,30cm
set termopt enhanced

set xlabel font ",30"
set ylabel font ",30"
set xtics font ",30"
set ytics font ",30"
set title font ",30"

set xlabel "{/Symbol D}Q"
set ylabel "T_{COM}"
set title "P=0.1, T_{imp}=0.1"

set xtics offset 0,-1
set xlabel offset 0,-3
set ylabel offset -3,0

set bmargin 10
set lmargin 20
set tmargin 10
set rmargin 20

set output "measurevscalc.pdf"
plot "gastemp_plot.dat" smooth unique w l lt rgb "blue" lw 4 title "Calculated",\
    "plotfile_clean.dat" using 1:2 smooth unique w l lt rgb "red" lw 4 title "Measured"
unset output
