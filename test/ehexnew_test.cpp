#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include "../include/functions.hpp"
#include "../include/globals.hpp"
#include "../include/classes.hpp"


int main(int argc, char** argv)
{
    setValues(0.2,0.01,0.2,4,0.1);

    double* rCM = new double[3];
    double* rCMtemp = new double[3];
    double* rCMstart = new double[3];
    double* vCM = new double[3];
    double* pos = new double[3];
    FILE* positions = fopen("ehexnew_test.xyz","w");
    FILE* tempout = fopen("ehextemp.dat","w");

    for(int i=0;i<3;i++)
        pos[i] = 0.;
    int run;

    Particle *cube = new Particle[N];
    std::list<Particle*> gas;
    InitPositions(cube);
    calcCM(cube,rCMStart,vCM);
    calcCM(cube,rCM,vCM);
    InitVelocities(cube);
    ComputeAccelerations(cube);

    for(run=0;run<40000;run++)
    {
            if(run%500==0)
                printf("(INIT) Zeitschritt %d\n",run);
            VelocityVerlet(cube,1,positions);
            if(run%100 == 0 && run < 20000)
                rescaleVelocities(cube);
            calcTemp(cube,tempout); 
    }
    for(run = 0;run<20000;run++)
    {
        if(run%500==0)
        {
            printf("(MEASURE) Zeitschritt %d - Number of Gas particles: %d\n",run,gas.size());
        }
        eHEXNew(cube);
        //BarostatNew(cube,gas);
        //harmonicTrap(rCM,vCM,rCMStart,cube);
        //if(run%400==0)
            //GenerateOutput(cube,gas,run+10000);
        //if(run%100==0)
            calcTemp(cube,tempout); 
    }
	return 0;
}
