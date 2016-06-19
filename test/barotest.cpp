#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include "../include/globals.hpp"
#include "../include/classes.hpp"
#include "../include/functions.hpp"

int main(int argc, char** argv)
{
    system("rm output/combined/*");
    FILE* tempdata = fopen("tempdata.dat","w");
    int run=0;
    gsl_rng_set(r,98);
    Particle *cube = new Particle[N];
    std::list<Particle*> gas;
    double* rCM = new double[3];
    double* vCM = new double[3];
    InitPositions(cube);
    calcCM(cube,rCMStart,vCM);
    calcCM(cube,rCM,vCM);
    InitVelocities(cube);
    ComputeAccelerations(cube);

    for(run=0;run<10000;run++)
    {
            if(run%200==0)
                printf("Zeitschritt %d\n",run);
            VelocityVerlet(cube,0,NULL);
            if(run%100 == 0)
                calcTemp(cube,tempdata); 
    }

    ComputeSoftSphere(gas,cube);

    for(run = 0;run<20000;run++)
    {
        if(run%100==0)
            printf("Zeitschritt %d\n",run);
        VelocityVerlet(cube,0,NULL);
        Barostat(cube,gas);
        harmonicTrap(rCM,vCM,rCMStart,cube);
        if(run%10==0)
            GenerateOutput(cube,gas,run);
        if(run%10==0)
            calcTemp(cube,tempdata); 
    }
    return 0;
}
