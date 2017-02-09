#!/bin/bash
#folders="$(ls -d p*)"
folders="$(ls -d *at0.1*)"

for i in $folders
do
    a="$(cat $i'/temperature_internal.dat' | wc -l)"
    if [ "$a" -ne 900 ] 
    then
        echo -e $i" (" $a ")"  
    fi
done
