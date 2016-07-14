#!/bin/bash
files="../*.dat"
#for f in $files
#do
    #echo $(basename ${f%.*})
#done
for f in $files
do
    filename=$(basename ${f%.*})
    gnuplot -e "\
        set term png; \
        set output '$filename.png'; \
        plot '$f' using 0:12 w l "
done

