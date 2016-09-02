#include "../include/globals.hpp"

//global variables
const gsl_rng_type *T = gsl_rng_taus;
gsl_rng *r = gsl_rng_alloc(T);
gsl_rng *r1 = gsl_rng_alloc(T);

int g_ID=1;
const unsigned int N = 864;

double* rCMStart = new double[3];
double rho = 1.1;
double Temp = 0.2;
double AmbientTemp = 0.9;
double P= 4.0;
double L = pow(N/rho,1.0/3);
double rCutOff = 2.5;
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
