set termopt enhanced;
set xlabel "Timestep";
set ylabel "T*";
set y2label "E_{tot}";
set ytics nomirror;
set y2tics;
set xrange [0:40000];
set x2range [0:40000];
set tics font ",15";
set xlabel font ",20";
set ylabel font ",20";
set y2label font ",20";
set key font ",15";
set key spacing 5;
plot "verlet_temp_1.dat" w l axes x1y1 title "Temperature" lw 2, "verlet_energy_1.dat" w l axes x2y2 title "Total Energy" lt rgb "blue" lw 2;

