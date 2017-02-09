#include "../include/globals.hpp"

//global variables
const gsl_rng_type *T = gsl_rng_taus;
gsl_rng *r = gsl_rng_alloc(T);
gsl_rng *r1 = gsl_rng_alloc(T);

gsl_histogram *gas_in = gsl_histogram_alloc(1000);
gsl_histogram *gas_out = gsl_histogram_alloc(1000);
gsl_histogram *gas_real_in = gsl_histogram_alloc(10000);

gsl_histogram2d *positionsxy = gsl_histogram2d_alloc(100,100);
gsl_histogram2d *positionsxz = gsl_histogram2d_alloc(100,100);
gsl_histogram2d *positionsyz = gsl_histogram2d_alloc(100,100);

int g_ID=1;
//const unsigned int N = 864;
//const unsigned int N = 500;
//const unsigned int N = 108;
//const unsigned int N = 256;
const unsigned int N = 32;

double* rCMStart = new double[3];
double* center = new double[3];
double rho = 0.9;
double Temp = 0.2;
double AmbientTemp = 0.9;
double P= 4.0;
double L = pow(N/rho,1.0/3);
double rCutOff = 2.5;
double rCut2 = rCutOff*rCutOff;


double dt = 0.01;
double eps = 0.5; //margin for boundary conditions
double dQ=0.050;
char* input;
//ACHTUNG: FG war bis jetzt nicht gesetzt!!! moegliche Aenderungen hier beginnen nachzuvollziehen
double FG = dQ/dt;

double*** Forces = new double**[N];
double*** Distances = new double**[N];

std::vector<double> virial;

bool PBC_FLAG = false;
bool EHEX_FLAG = false;


//enum Cubeface { X0, XL, Y0, YL, Z0, ZL };
