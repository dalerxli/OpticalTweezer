#!/bin/bash
cd ~/codes/OpticalTweezer/output/runs
folders="$(ls -d */ | tail -n 9 | head -n 7)"
#folders="$(ls -d ~/codes/OpticalTweezer/output/runs/*/ | tail -n 9 | head -n 7)"
for i in $folders
do
    echo $i;
    cat $i"temperature_internal.dat" | wc -l;
done
