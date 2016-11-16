#!/bin/bash

g++  cell_list.cpp ../src/functions.cpp ../src/globals.cpp ../src/classes.cpp -lgsl -lgslcblas  -o  cell_list.o
