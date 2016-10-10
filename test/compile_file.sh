#!/bin/bash

g++ comtest.cpp ../src/!(main).cpp -lgsl -lgslcblas -o comtest.o
