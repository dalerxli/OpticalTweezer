#!/bin/bash

files=*.dat
directory=changed/
#files="$(ls | sort --version-sort -f | head -n 1000)"
#rm grepout.dat
if [ ! -d $directory ]
then
    echo "Creating directory" $directory
    mkdir $directory
else
    numfiles="$(find changed/ -type f -name '*.dat' | wc -l)"
    if [ "$numfiles" -ne 0 ]
    then
        #rm changed/*.dat
        subfiles=changed/*.dat
        for i in $subfiles
        do
            echo "removing" $i
            rm $i
        done
    fi
fi
#files="$(ls | sort --version-sort -f | tail -n $((numfiles-10)) | head -n 10)"
#echo $files
count=0
for i in $files
do
    if [ "$((count % 1000))" -eq 0 ]
    then
        echo $i
    fi
    #echo $i >> grepout.dat
    #echo "$(cat $i | grep "GasChange")" >> grepout.dat
    $(cat $i | grep -q 'GasChange')
    if [ "$?" -eq 0 ] 
    then
        #echo "File" $i "does not contain GasChange"
        cp $i changed/
    fi
    #count=$((count+1))
done
#yes="$(cat 865.dat | grep 'GasChange')"
#echo $yes

