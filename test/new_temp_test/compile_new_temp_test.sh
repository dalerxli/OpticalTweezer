#!/bin/bash
rm combined/*
g++  new_temp_test.cpp ../../src/functions.cpp ../../src/globals.cpp ../../src/classes.cpp -lgsl -lgslcblas -o  new_temp_test.o
