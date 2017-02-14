set term pdfcairo size 40cm,30cm
set termopt enhanced
set xlabel font ",50"
set ylabel font ",50"
set tics font ",50"
set title font ",50"
set key off

set xlabel offset 0,-4
set ylabel offset -5,0
set xtics offset 0,-1
set title offset 0,1


set title "P=0.5"
set xlabel "{/Symbol D}Q"
set ylabel "T_{INT}"

set yrange [0:1.4]
set xrange [0:*]

set bmargin 10
set lmargin 20
set tmargin 10
set rmargin 20

set output "../doc/images/p05_int.pdf"
plot "p0.5_at0.05/plotfile_clean.dat" using 1:4 smooth unique with linespoints pt 7  lw 4 lt rgb "red" title "T_{imp} = 0.05", \
    "p0.5_at0.1/plotfile_clean.dat" using 1:4 smooth unique with linespoints pt 7  lw 4 lt rgb "blue" title "T_{imp} = 0.1"
unset output
