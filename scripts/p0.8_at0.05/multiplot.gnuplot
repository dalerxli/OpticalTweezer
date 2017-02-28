set term pdfcairo size 40cm,30cm
set termopt enhanced
set title "P=0.8    T_{imp}=0.05"
set title font ",30"
set xlabel "Time"
set ylabel "Temperature"
set xlabel font ",30"
set xlabel offset 0,-5
set ylabel font ",30"
set ylabel offset -5,0
set xtics font ",30"
set xtics offset 0,-2
set ytics font ",30"
set key font ",30"
set key spacing 5
set bmargin 10
set lmargin 20
set tmargin 10
set rmargin 20
#set output "multiplot.pdf"
#set multiplot layout 1,3 title "P=0.8    T_{imp}=0.05"
set output "t_int.pdf"
plot "p0.8_at0.05_170206_1132_temperature_internal.dat" using ($0*0.01*100):1 w l lw 3 title "T_{int}",0.048 lw 6 lt rgb "blue" title "<T_{int}>";
unset output
set output "gastemp.pdf"
plot "p0.8_at0.05_170206_1132_gasTempData.dat" using ($0*0.01*100):2 w l lw 3 title "T_{em}",0.0494 lw 6 lt rgb "blue" title "<T_{em}>";
unset output
set output "vcomtemp.pdf"
set xlabel "Velocity"
set ylabel "Count"
plot "p0.8_at0.05_170206_1132_histogram_mass_vx.dat" with boxes lw 4 lt rgb "red" title "v_x",sqrt(32/(2*pi*0.051))*exp(-x**2*32/(2*0.051)) lw 6 lt rgb "blue" title "p(v_x)";
unset output
set output "vcomsqd.pdf"
set xlabel "Time"
set ylabel "Temperature"
set ylabel offset -5,0
plot "p0.8_at0.05_170206_1132_vCOMData.dat" using ($0*0.01*100):(($1*$1+$2*$2+$3*$3)*0.5*32) with lines lw 4 lt rgb "red" title "T_{COM}",0.051 lw 6 lt rgb "blue" title "<T_{COM}>";
unset output
