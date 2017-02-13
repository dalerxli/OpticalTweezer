#!/bin/bash
rm [0-9]*

cd ~/GitHub/OpticalTweezer/test/output/runs/
#folders="$(ls -d */ | tail -n 9 | head -n 7)"
folders=$1

for i in $folders
do
    #cp $i"vCOMData.dat" ../../scripts/data/"vCOMData_"${i%/}".dat";
    #cp $i"parameters.txt" ../../scripts/data/"parameters_"${i%/}".txt";
    #cp $i"gasTempData.dat" ../../scripts/data/"gasTempData"${i%/}".dat";
    cp $i"/vCOMData.dat" ../../scripts/data3/${i%/}"_vCOMData.dat";
    cp $i"/parameters.txt" ../../scripts/data3/${i%/}"_parameters.txt";
    cp $i"/gasTempData.dat" ../../scripts/data3/${i%/}"_gasTempData.dat";

done
