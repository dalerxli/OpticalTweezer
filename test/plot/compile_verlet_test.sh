#!/bin/bash
sources="$(cat sources.txt)"
g++  verlet_test.cpp  $sources -lgsl -lgslcblas -o  verlet_test.o
