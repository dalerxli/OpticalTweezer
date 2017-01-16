#!/bin/bash
sources="$(cat sources.txt)"
g++  baro_output_equib.cpp  $sources -lgsl -lgslcblas -o  baro_output_equib.o
