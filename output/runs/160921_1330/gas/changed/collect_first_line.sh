#!/bin/bash

#files="$(find ./ -name '*.dat')"
files=*.dat
numfiles="$(find ./ -name '*.dat' | wc -l)"
script="/home/motz/codes/python/calc_percent.py"
outfile=first_lines.txt
if [ -e $outfile ] 
then
    rm $outfile
fi
count=0
for i in $files
do
    if [ "$((count % 1000))" -eq 0 ]
    then
        number="$(echo -e $count"\n"$numfiles"\n" | python $script)"
        echo $count "of" $numfiles "("$number"%)"
    fi
    cat $i | head -n 1 >> $outfile
    count=$((count+1))
done
