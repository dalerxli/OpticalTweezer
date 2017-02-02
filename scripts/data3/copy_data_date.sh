#!/bin/bash
datetime=$1

cp ../../test/output/runs/$datetime"/vCOMData.dat" ./${datetime%/}"_vCOMData.dat";
cp ../../test/output/runs/$datetime"/parameters.txt" ./${datetime%/}"_parameters.txt";
cp ../../test/output/runs/$datetime"/gasTempData.dat" ./${datetime%/}"_gasTempData.dat";
