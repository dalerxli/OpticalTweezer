#!/bin/bash

datetime="$(ls | grep -P -o "[0-9]{6}_[0-9]{4}" | sort -u)"
for i in $datetime
do
    python vCOMhistogram.py $i
    ./create_gnuplot.sh $i
done


