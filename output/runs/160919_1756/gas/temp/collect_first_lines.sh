#!/bin/bash

rm values.dat
files=*.dat
counter=0
for i in $files
do
    if [ $((counter % 2)) -eq 0 ] 
    then
        echo $counter
    fi
    head -n 1 $i | awk '{print $4}' >> values.dat
    counter=$((counter+1))
done
