#!/bin/bash
sources="$(cat sources.txt)"
g++  comhist.cpp  $sources -lgsl -lgslcblas -o  comhist.o
