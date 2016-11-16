#!/bin/bash

g++  temperature_test.cpp ../../src/functions.cpp ../../src/globals.cpp ../../src/classes.cpp -lgsl -lgslcblas  -o  temperature_test.o
