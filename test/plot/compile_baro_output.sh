#!/bin/bash
sources="$(cat sources.txt)"
g++  baro_output.cpp  $sources -lgsl -lgslcblas -o  baro_output.o
