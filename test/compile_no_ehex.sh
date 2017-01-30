#!/bin/bash
sources="$(cat sources.txt)"
g++  no_ehex_test.cpp  $sources -lgsl -lgslcblas -o bin/no_ehex_test.o
