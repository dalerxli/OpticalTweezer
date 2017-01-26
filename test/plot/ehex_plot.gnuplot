set term pdfcairo size 40cm,30cm
set termopt enhanced
set key at 6500,3
set key spacing 6
set key nobox
set xlabel font ",30"
set ylabel font ",30"
set xlabel "Time"
set ylabel "T*"
set xlabel offset 0,-7
set ylabel offset -7,0
set bmargin 15
set lmargin 20
set tmargin 10
set rmargin 20
set key font ",30"
set tics font ",30"
set xtics offset 0,-1
set palette rgb 33,13,10;
unset colorbox
set output "ehex_temp_dq.pdf";
plot for [i=0:9] "ehex_temp_00".i."_1.dat" w l title "{/ Symbol D}Q = 0.0".i lw 4 lt palette cb i;
#plot for [i=0:9] "ehex_temp_00".i."_1.dat" using ($0*0.0217):1 w l title "{/ Symbol D}Q = 0.0".i lw 4 lt palette cb i;
unset output;
#set output "ehex_temp_dq.pdf";
#plot "ehex_temp_000_1.dat" w l lw 2 lt 1 lc rgb "#00ff00"
#unset output;


