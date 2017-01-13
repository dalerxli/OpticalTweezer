set term pdfcairo size 40cm,30cm
set key left
set xlabel font ",30"
set ylabel font ",30"
set xlabel "Timestep"
set ylabel "T*"
set xlabel offset 0,-7
set ylabel offset -7,0
set bmargin 15
set lmargin 20
set tmargin 10
set rmargin 20
set key font ",30"
set key spacing 1.5
set tics font ",30"
set palette rgb 33,13,10;
unset colorbox
set output "ehex_temp_dq.pdf";
plot for [i=0:9] "ehex_temp_00".i."_1.dat" w l title "q = 0.0".i lw 2 lt palette cb i;
unset output;
#set output "ehex_temp_dq.pdf";
#plot "ehex_temp_000_1.dat" w l lw 2 lt 1 lc rgb "#00ff00"
#unset output;


