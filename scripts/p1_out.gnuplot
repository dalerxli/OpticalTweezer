set term pdfcairo size 40cm,30cm
set termopt enhanced
set xlabel font ",30"
set ylabel font ",30"
set key font ",30"
set tics font ",30"
set title font ",30"

set xlabel offset 0,-3
set ylabel offset -5,0
set xtics offset 0,-1
set title offset 0,1

set key spacing 5

set title "P=1.0"
set xlabel "{/Symbol D}Q"
set ylabel "T_{em}"
set yrange [0:0.13]
set xrange [0:0.5]

set bmargin 10
set lmargin 20
set tmargin 10
set rmargin 20

set output "p1_out.pdf"
plot "p1_at0.05/plotfile_clean.dat" using 1:8 smooth unique lw 4 lt rgb "red" title "T_{em} = 0.05", \
    "p1_at0.1/plotfile_clean.dat" using 1:8 smooth unique lw 4 lt rgb "blue" title "T_{em} = 0.1"
unset output

