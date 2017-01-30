#!/bin/bash
rm combined/*
sources="$(cat sources.txt)"
g++  barotemp_new.cpp  $sources -lgsl -lgslcblas -o  barotemp_new.o
