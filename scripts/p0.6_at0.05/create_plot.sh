#!/bin/bash
if [ -e "plotfile.dat" ]
then 
    rm plotfile.dat
fi
if [ -e "plotfile_clean.dat" ]
then 
    rm plotfile_clean.dat
fi
datetimes="$(ls * | grep -P -o "p[0-9].?[0-9]{0,2}_at[0-9].[0-9]{1,2}_[0-9]{6}_[0-9]{4}" | sort -u)"

for i in $datetimes
do
    q="$(cat "$i"_parameters.txt | grep "dQ" | cut -f3 -d " ")"
    tcom="$(cat "$i"_statistics_mass.dat | grep "COM" | cut -f2 -d " ")"
    terr="$(cat "$i"_statistics_mass.dat | grep "COM" | cut -f4 -d " ")"
    ti="$(cat "$i"_statistics_mass.dat | grep "INT" | cut -f2 -d " ")"
    interr="$(cat "$i"_statistics_mass.dat | grep "INT" | cut -f4 -d " ")"
    gin="$(cat "$i"_statistics_mass.dat | grep "In" | cut -f2 -d " ")"
    ginerr="$(cat "$i"_statistics_mass.dat | grep "In" | cut -f4 -d " ")"
    gout="$(cat "$i"_statistics_mass.dat | grep "Out" | cut -f2 -d " ")"
    gouterr="$(cat "$i"_statistics_mass.dat | grep "Out" | cut -f4 -d " ")"
    echo -e $q"\t\t"$tcom"\t\t"$terr"\t\t"$ti"\t\t"$interr"\t\t"$gin"\t\t"$ginerr"\t\t"$gout"\t\t"$gouterr>> plotfile.dat
done

cp plotfile.dat plotfile_clean.dat
