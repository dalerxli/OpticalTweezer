#!/bin/bash
dir="$(pwd)"
target="$(basename $dir)"
#target=$1
folders="$(ls -d ../../test/output/runs/$target*)"
if [ -e folderlist.txt ]
then 
    rm folderlist.txt
fi
for i in $folders
do
    basename $i >> folderlist.txt
done

names="$(cat folderlist.txt)"

for i in $names
do
    ./copy_data_date.sh $i
done

for i in $names
do
    python vCOMhistogramMass.py $i
done

./create_png.sh

gnuplot -e "\
    files=system('ls *.gnuplot');
    do for [ i in files ] {
        load i
    }"

./create_plot.sh
