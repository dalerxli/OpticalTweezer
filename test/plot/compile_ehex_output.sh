#!/bin/bash
sources="$(cat sources.txt)"
g++  ehex_output.cpp  $sources -lgsl -lgslcblas -o  ehex_output.o
