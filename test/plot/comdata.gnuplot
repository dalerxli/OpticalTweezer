# plot data for vCOM 
# data used: 170111_1241
folder="~/codes/OpticalTweezer/output/runs/170111_1241/"
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
set ylabel "r_{COM}"
set xlabel "Timestep"
set output "rCOM_nosquare.pdf"
plot folder."comdata.dat" using 1 w l lw 5 title "r_x", "" using 2 w l lw 5title "r_y", "" using 3 w l lw 5 title "r_z"
unset output

