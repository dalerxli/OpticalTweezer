#!/bin/bash

files="*.dat"
numfiles="$(ls -1 | wc -l)"
mean=0
for i in $files
do
    number="$(cat $i | wc -l)"
    mean=$((mean+number))
done
mean=$((mean/numfiles))
echo $mean
