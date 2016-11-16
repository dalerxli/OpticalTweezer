#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include "../include/globals.hpp"
#include "../include/classes.hpp"
#include "../include/functions.hpp"


struct cell{
    double xStart;
    double yStart;
    double zStart;
};
int main(int argc, char** argv)
{
    int a;
    int run=0;
    gsl_rng_set(r,98);
    gsl_histogram_set_ranges_uniform(gas_in,-1,5);
    gsl_histogram_set_ranges_uniform(gas_out,-1,5);
    gsl_histogram_set_ranges_uniform(gas_real_in,-1,5);

    Particle *cube = new Particle[N];
    std::list<Particle*> gas;
    std::list<Particle*> gasHistory;
    double* rCM = new double[3];
    double* rCMtemp = new double[3];
    double* vCM = new double[3];
    double energy = 0;

    InitPositions(cube);
    calcCM(cube,rCMStart,vCM);
    calcCM(cube,rCM,vCM);
    InitVelocities(cube);
    ComputeAccelerations(cube);
    /*
     *for(run=0;run<10000;run++)
     *{
     *        if(run%500==0)
     *            printf("(INIT) Zeitschritt %d\n",run);
     *        VelocityVerlet(cube,0,NULL);
     *        if(run%100 == 0)
     *            rescaleVelocities(cube);
     *        if(run%100 == 0)
     *            GenerateOutput(cube,gas,run);
     *}
     */

    gsl_histogram *cells = gsl_histogram_alloc(100);
    gsl_histogram_set_ranges_uniform(cells,0,50);
    double lcx = floor(L/rCutOff);
    double rcx = L/lcx;
    double c = 0;
    double ci[3] = {0,0,0};
    for(int i=0;i<N;i++)
    {
        for(int k=0;k<3;k++)
            ci[k] = floor(cube[i].r[k]/rcx);
        c = ci[0]*lcx*lcx+ci[1]*lcx+ci[2];
        std::cout << c << std::endl;
        gsl_histogram_increment(cells,c);
    } 

    FILE* celldata = fopen("cell_list_test.dat","w");
    gsl_histogram_fprintf(celldata,cells,"%f","%f");
    fclose(celldata);



	return 0;
}
