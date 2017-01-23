#!/bin/bash
sources="$(cat sources.txt)"
g++  newalgorithm.cpp $sources -lgsl -lgslcblas -o  newalgorithm.o
