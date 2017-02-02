#!/bin/bash
datetime=$1
#echo -e "This is "$datetime
echo -e 'Tin=system("./get_number.sh GasIn '$datetime'_statistics.dat")+0' > $datetime"_plot.gnuplot"
echo -e 'Tout=system("./get_number.sh GasOut '$datetime'_statistics.dat")+0' >> $datetime"_plot.gnuplot"
echo -e 'Tcm=system("./get_number.sh T_COM '$datetime'_statistics.dat")+0' >> $datetime"_plot.gnuplot"
echo -e 'plot "'$datetime'_histogram_vx.dat" with boxes, \
    sqrt(1./(2.*pi*Tcm))*exp(-x**2/(2*Tcm)), \
    sqrt(0.1/(2.*pi*Tin))*exp(-x**2*0.1/(2*Tin)), \
    sqrt(0.1/(2.*pi*Tout))*exp(-x**2*0.1/(2*Tout))' >> $datetime"_plot.gnuplot"
    #sqrt(1./(2.*pi*Tin))*exp(-x**2/(2*Tin)), \
    #sqrt(1./(2.*pi*Tout))*exp(-x**2/(2*Tout))' >> $datetime"_plot.gnuplot"


