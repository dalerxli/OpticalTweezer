#ifndef _INCLUDE_GLOBALS_HPP_
#define _INCLUDE_GLOBALS_HPP_

#include <iostream>
#include <cstdio>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <list>
#include <iterator>
#include <vector>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cstring>
#include <fstream>

#define _TEST printf("works!\n") 

//global variables
extern gsl_rng *r;
extern gsl_rng *r1;
extern const gsl_rng_type *T;

extern int g_ID;
extern const unsigned int N;

extern double FG;
extern double* rCMStart;
extern const double eps;
extern const double dQ;
extern const double rho;
extern const double L;
extern const double rCutOff;
extern const double dt;

extern bool EHEX_FLAG;
extern bool PBC_FLAG;

enum Cubeface { X0, XL, Y0, YL, Z0, ZL };

#endif
