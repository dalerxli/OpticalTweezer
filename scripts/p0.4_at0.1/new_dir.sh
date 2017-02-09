#!/bin/bash

mkdir $1
cp *.sh $1/
cp *.py $1/
cd $1
./init.sh
