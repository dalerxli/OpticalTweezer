#!/bin/bash

g++  outputs.cpp ../src/functions.cpp ../src/globals.cpp ../src/classes.cpp -lgsl -lgslcblas -o  outputs.o
