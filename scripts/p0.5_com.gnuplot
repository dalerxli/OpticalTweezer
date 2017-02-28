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

set yrange [0:0.5]
set xrange [0:*]


set title "P=0.5"
set xlabel "{/Symbol D}Q"
set ylabel "T_{COM}"

set bmargin 10
set lmargin 30
set tmargin 10
set rmargin 30


set output "../doc/images/p05_com.pdf"
plot "p0.5_at0.05/plotfile_clean.dat" using 1:2 smooth unique with linespoints pt 20 lw 25 lt rgb "red" title "T_{imp} = 0.05", \
    "p0.5_at0.1/plotfile_clean.dat" using 1:2 smooth unique with linespoints pt 20 lw 25 lt rgb "blue" title "T_{imp} = 0.1"
unset output
