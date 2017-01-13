set term pdfcairo size 50cm,40cm
set termopt enhanced;
set xlabel "Timestep";
set ylabel "T*";
set y2label "E_{tot}";
set ytics nomirror;
set y2tics;
set xrange [0:30000];
set x2range [0:30000];
set yrange [0:0.4];
set tics font ",30";
set xlabel font ",30";
set ylabel font ",30";
set y2label font ",30";
set key font ",30";
set key spacing 2;
set lmargin 20;
set tmargin 10;
set rmargin 20;
set bmargin 20;
set key spacing 2
set key font ",30"
set xtics offset 0,-2
set ylabel offset -10,0
set y2label offset 10,0
set rmargin 30
set lmargin 30
set bmargin 25
set xlabel offset 0,-5
set tics font ",30"
set output "verlet_output_te.pdf";plot "verlet_temp_1.dat" w l axes x1y1 title "Temperature" lw 2 lt rgb "red", "verlet_energy_1.dat" w l axes x2y2 title "Total Energy" lt rgb "blue" lw 2;unset output;
#plot "verlet_temp_1.dat" w l axes x1y1 title "Temperature" lw 2 lt rgb "red", "verlet_energy_1.dat" w l axes x2y2 title "Total Energy" lt rgb "blue" lw 2;

### commands for fine tuning:
#   set key spacing 2
#   set key font ",30"
#   set xtics offset 0,-2
#   set ylabel offset -10,0
#   set y2label offset 10,0
#   set rmargin 30
#   set lmargin 30
#   set bmargin 25
#   set xlabel offset 0,-5
#   set tics font ",30"
#   set term pdfcairo size 50cm,40cm
#   set output "verlet_output_te.pdf";plot "verlet_temp_1.dat" w l axes x1y1 title "Temperature" lw 2 lt rgb "red", "verlet_energy_1.dat" w l axes x2y2 title "Total Energy" lt rgb "blue" lw 2;unset output;
