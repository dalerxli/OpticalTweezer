set term pdfcairo size 40cm,30cm
set termopt enhanced

set xlabel font ",30"
set ylabel font ",30"
set xtics font ",30"
set ytics font ",30"
set title font ",30"
set key font ",30"
set key spacing 10
set key left at 0.04,0.54

set xlabel "{/Symbol D}Q"
set ylabel "T_{COM}"
set title "P=0.8, T_{imp}=0.05"

set xtics offset 0,-1
set xlabel offset 0,-3
set ylabel offset -3,0

set bmargin 10
set lmargin 20
set tmargin 10
set rmargin 20

set output "measurevscalc.pdf"
plot "gastemp_plot.dat" smooth unique with lines lt rgb "blue" lw 10 title "Calculated",\
    "plotfile_clean.dat" using 1:2 with linespoints lt rgb "red" pt 20 lw 10 title "Measured"
unset output
