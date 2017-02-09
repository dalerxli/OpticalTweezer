#!/bin/bash
var=$1
file=$2
grep -o -P -e $var":.*" $file | cut -f2- -d" "
#grep -o -P -e "GasIn:.*" 170123_2033_statistics.dat | cut -f2- -d:
# gets number after GasIn in file 170123_2033_statistics.dat
