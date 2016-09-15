#!/bin/bash

#i=12
#mkdir $i
#cp ../$i* $i/
for i in {10..14}
do
    mkdir $i
    cp ../$i* $i/
    sed -i '/test/d' $i/*
    sed -i '/865/d' $i/*
    sed -i '/Glass/d' $i/*
done
