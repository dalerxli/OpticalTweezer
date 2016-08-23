#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include "../include/classes.hpp"
#include "../include/functions.hpp"
#include "../include/globals.hpp"

int main(int argc, char** argv)
{
    gsl_rng_set(r,98);
    FILE* tempout = fopen("temperature_internal.dat","w");
    FILE* comData = fopen("comdata.dat","w");
    Particle *cube = new Particle[N];
    std::list<Particle*> gas;
    double* rCM = new double[3];
    double* rCMtemp = new double[3];
    double* vCM = new double[3];
    double energy = 0;
    int run;
    InitPositions(cube);
    calcCM(cube,rCMStart,vCM);
    calcCM(cube,rCM,vCM);
    InitVelocities(cube);
    ComputeAccelerations(cube);
    for(run=0;run<10000;run++)
    {
            if(run%500==0)
                printf("(INIT) Zeitschritt %d\n",run);
            VelocityVerlet(cube,0,NULL);
            if(run%100 == 0)
            {
                rescaleVelocities(cube);
                calcTemp(cube,tempout); 
            }
            if(run%400 == 0)
                GenerateOutput(cube,gas,run);
    }
    ComputeSoftSphere(gas,cube);

    for(run = 0;run<40000;run++)
    {
        if(run%500==0)
            printf("(MEASURE)Zeitschritt %d\n",run);
        eHEX(cube);
        BarostatNew(cube,gas);
        //std::cout << "works!" << std::endl;
        harmonicTrap(rCM,vCM,rCMStart,cube);
        //std::cout << "works!" << std::endl;
        if(run%400==0)
        {
            //trackParticle(cube,gas,1400,particleTracker);
            GenerateOutput(cube,gas,run+10000);
            //std::cout << "works!" << std::endl;
            calcCM(cube,rCMtemp,comData);
        }
        if(run%100==0)
            calcTemp(cube,tempout); 
    }
    chdir("../../../");
    delete [] cube;
    gas.clear(); 
    g_ID = 1;
	EHEX_FLAG = false;

    fclose(tempout);
    fclose(comData);
    return 0;
}
