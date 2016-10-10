#!/bin/bash

files="$(ls | sort --version-sort -f | head -n 1000)"


for i in $files
do
    grep -q -l "GasChange" $i
    if [ "$?" -eq 0 ] 
    then
        cp $i changed2/
    fi
done 
