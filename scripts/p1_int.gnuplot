set term pdfcairo size 40cm,30cm
set termopt enhanced
set xlabel font ",70"
set ylabel font ",70"
set tics font ",70"
set title font ",70"
set key off

set xlabel offset 0,-5
set ylabel offset -14,0
set xtics offset 0,-2
set title offset 0,1


set title "P=1.0"
set xlabel "{/Symbol D}Q"
set ylabel "T_{INT}"

set yrange [0:1.4]
set xrange [0:*]

set bmargin 10
set lmargin 30
set tmargin 10
set rmargin 30


set output "../doc/images/p1_int.pdf"
plot "p1_at0.05/plotfile_clean.dat" using 1:4 smooth unique with linespoints pt 20 lw 25 lt rgb "red" title "T_{imp} = 0.05", \
    "p1_at0.1/plotfile_clean.dat" using 1:4 smooth unique with linespoints pt 20 lw 25 lt rgb "blue" title "T_{imp} = 0.1"
unset output
