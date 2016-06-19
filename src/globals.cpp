#include "../include/globals.hpp"

//global variables
const gsl_rng_type *T = gsl_rng_taus;
gsl_rng *r = gsl_rng_alloc(T);
gsl_rng *r1 = gsl_rng_alloc(T);

int g_ID=1;
const unsigned int N = 864;

double FG;
double* rCMStart = new double[3];
const double rho = 1.1;
const double L = pow(N/rho,1.0/3);
const double rCutOff = 2.5;
const double dt = 0.01;
const double eps = 0.5;
const double dQ=0.010;

bool PBC_FLAG = false;
bool EHEX_FLAG = false;

//enum Cubeface { X0, XL, Y0, YL, Z0, ZL };
