#!/bin/bash
sources="$(cat sources.txt)"
g++  whereisthecom.cpp $sources -lgsl -lgslcblas -o  whereisthecom.o
