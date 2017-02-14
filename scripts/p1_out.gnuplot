set term pdfcairo size 40cm,30cm
set termopt enhanced
set xlabel font ",50"
set ylabel font ",50"
set tics font ",50"
set title font ",50"
set key off

set xlabel offset 0,-4
set ylabel offset -5,2
set xtics offset 0,-1
set title offset 0,1


set title "P=1.0"
set xlabel "{/Symbol D}Q"
set ylabel "T_{em}"
set yrange [0:0.13]
set xrange [0:*]

set bmargin 10
set lmargin 20
set tmargin 10
set rmargin 20

set output "../doc/images/p1_out.pdf"
plot "p1_at0.05/plotfile_clean.dat" using 1:8 smooth unique with linespoints lw 4 pt 7 lt rgb "red" title "T_{em} = 0.05", \
    "p1_at0.1/plotfile_clean.dat" using 1:8 smooth unique with linespoints lw 4 pt 7 lt rgb "blue" title "T_{em} = 0.1", \
    #"p1_at0.05/plotfile_clean.dat" using 1:8 with points pointtype 5 notitle, \
    #"p1_at0.1/plotfile_clean.dat" using 1:8 with points pointtype 5 notitle
unset output

