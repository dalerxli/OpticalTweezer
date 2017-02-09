#!/bin/bash
if [ -e "plotfile.dat" ]
then 
    rm plotfile.dat
fi
datetimes="$(ls * | grep -P -o "p[0-9].?[0-9]{0,2}_at[0-9].[0-9]{1,2}_[0-9]{6}_[0-9]{4}" | sort -u)"

for i in $datetimes
do
    q="$(cat "$i"_parameters.txt | grep "dQ" | cut -f3 -d " ")"
    tcom="$(cat "$i"_statistics_mass.dat | grep "COM" | cut -f2 -d " ")"
    ti="$(cat "$i"_statistics_mass.dat | grep "INT" | cut -f2 -d " ")"
    gin="$(cat "$i"_statistics_mass.dat | grep "In" | cut -f2 -d " ")"
    gout="$(cat "$i"_statistics_mass.dat | grep "Out" | cut -f2 -d " ")"
    echo -e $q"\t"$tcom"\t"$ti"\t"$gin"\t"$gout >> plotfile.dat
done
