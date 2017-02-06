#!/bin/bash

datetimes="$(ls * | grep -P -o "p[0-9]{2}_[0-9]{6}_[0-9]{4}" | sort -u)"

for i in $datetimes
do
    plotfile=$i"_pngplot.gnuplot"
    echo -e "set term png \n" > $plotfile
    echo -e "b = system('cat " $i"_parameters.txt') \n" >> $plotfile
    echo -e "set title b\n" >> $plotfile
    echo -e "set output '"$i"_temperature_internal.png'\n" >> $plotfile
    echo -e "plot '"$i"_temperature_internal.dat' w l \n" >> $plotfile
    echo -e "unset output \n" >> $plotfile
done
