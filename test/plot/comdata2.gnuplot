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
set key spacing 1.5
set ylabel "r_{COM}^2"
set xlabel "Timestep"
set output "comdata.pdf"
plot folder."comdata.dat" using ($1*$1+$2*$2+$3*$3) w l lw 5 title "r_{COM}^2"
unset output

