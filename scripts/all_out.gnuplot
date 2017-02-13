set term pdfcairo size 40cm,30cm
set termopt enhanced
set xlabel font ",30"
set ylabel font ",30"
set key font ",30"
set tics font ",30"
set title font ",30"

set palette rgb 30,31,32
set cbrange [1:7]
unset colorbox


set xlabel offset 0,-3
set ylabel offset -5,0
set xtics offset 0,-1
set title offset 0,1

set key spacing 5
set key right at 0.5,0.129

set title "Emerging Gas Temperature"
set xlabel "{/Symbol D}Q"
set ylabel "T_{em}"

set yrange [0:0.13]
set xrange [0:0.5]

set bmargin 10
set lmargin 20
set tmargin 10
set rmargin 20

set output "all_out.pdf"
plot "p0.1_at0.05/plotfile_clean.dat" using 1:8 smooth unique lw 4 lt palette cb 1 title "P = 0.1, T_{em} = 0.05", \
    "p0.1_at0.1/plotfile_clean.dat" using 1:8 smooth unique lw 4 lt palette cb 2 title "P = 0.1, T_{em} = 0.1",\
    "p0.5_at0.05/plotfile_clean.dat" using 1:8 smooth unique lw 4 lt palette cb 3 title "P = 0.5, T_{em} = 0.1",\
    "p0.5_at0.1/plotfile_clean.dat" using 1:8 smooth unique lw 4 lt palette cb 4 title "P = 0.5, T_{em} = 0.1",\
    "p1_at0.05/plotfile_clean.dat" using 1:8 smooth unique lw 4 lt palette cb 5 title "P = 1.0, T_{em} = 0.1",\
    "p1_at0.1/plotfile_clean.dat" using 1:8 smooth unique lw 4 lt palette cb 6 title "P = 1.0, T_{em} = 0.1"
unset output

