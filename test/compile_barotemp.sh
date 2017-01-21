#!/bin/bash
rm combined/*
sources="$(cat sources.txt)"
g++  barotemp.cpp  $sources -lgsl -lgslcblas -o  barotemp.o
