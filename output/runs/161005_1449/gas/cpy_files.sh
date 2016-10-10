#!/bin/bash

files=*.dat
count=0
for i in $files
do
    cp $i copy/
    if [ "$((count%1000))" -eq 0 ]
    then
        echo $count
    fi
    count=$((count+1))
done
