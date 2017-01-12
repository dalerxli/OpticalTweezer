#!/bin/bash
sources="$(cat sources.txt)"
g++  verlet_output.cpp  $sources -lgsl -lgslcblas -o  verlet_output.o
