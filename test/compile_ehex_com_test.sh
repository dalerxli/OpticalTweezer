#!/bin/bash
sources="$(cat sources.txt)"
g++  ehex_com_test.cpp  $sources -lgsl -lgslcblas -o bin/ehex_com_test.o

