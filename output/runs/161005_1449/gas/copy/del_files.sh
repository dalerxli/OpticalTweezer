#!/bin/bash

files=*.dat
count=0
for i in $files
do
    if [ "$((count%1000))" -eq 0 ]
    then 
        echo $count
    fi
    grep -e "GasChange" -e "GasOut" -q $i 
    if [ "$?" -eq 0 ]
    then
        rm $i
    fi
    count=$((count+1))
done
