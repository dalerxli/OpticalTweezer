#!/bin/bash
sources="$(cat sources.txt)"
g++  inertia_test.cpp $sources -lgsl -lgslcblas -o  inertia_test.o
