#include "../include/globals.hpp"

//global variables
const gsl_rng_type *T = gsl_rng_taus;
gsl_rng *r = gsl_rng_alloc(T);
gsl_rng *r1 = gsl_rng_alloc(T);

int g_ID=1;
const unsigned int N = 864;

double* rCMStart = new double[3];
double rho = 1.1;
double Temp;
double L = pow(N/rho,1.0/3);
double rCutOff = 2.5;
double dt = 0.01;
double eps = 0.5;
double dQ=0.050;
char* input;
//ACHTUNG: FG war bis jetzt nicht gesetzt!!! moegliche Aenderungen hier beginnen nachzuvollziehen
double FG = dQ/dt;

bool PBC_FLAG = false;
bool EHEX_FLAG = false;

//enum Cubeface { X0, XL, Y0, YL, Z0, ZL };
