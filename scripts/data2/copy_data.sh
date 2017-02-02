#!/bin/bash
rm [0-9]*

cd ~/codes/OpticalTweezer/output/runs
folders="$(ls -d */ | tail -n 8)"

for i in $folders
do
    #cp $i"vCOMData.dat" ../../scripts/data/"vCOMData_"${i%/}".dat";
    #cp $i"parameters.txt" ../../scripts/data/"parameters_"${i%/}".txt";
    #cp $i"gasTempData.dat" ../../scripts/data/"gasTempData"${i%/}".dat";
    cp $i"vCOMData.dat" ../../scripts/data2/${i%/}"_vCOMData.dat";
    cp $i"parameters.txt" ../../scripts/data2/${i%/}"_parameters.txt";
    cp $i"gasTempData.dat" ../../scripts/data2/${i%/}"_gasTempData.dat";
done
