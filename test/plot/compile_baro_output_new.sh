#!/bin/bash
sources="$(cat sources.txt)"
g++  baro_output_new.cpp  $sources -lgsl -lgslcblas -o  baro_output_new.o
