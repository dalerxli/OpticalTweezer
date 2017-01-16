#!/bin/bash
sources="$(cat sources.txt)"
g++  baro_output_no_ehex.cpp  $sources -lgsl -lgslcblas -o  baro_output_no_ehex.o

