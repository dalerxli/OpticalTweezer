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
    double* rCM = new double[3];
    double* rCMtemp = new double[3];
    double* rCMstart = new double[3];
    double* vCM = new double[3];
    double* pos = new double[3];
    FILE* positions = fopen("harmonicTrapData.xyz","w");

    for(int i=0;i<3;i++)
        pos[i] = 0.;
    int run;

    Particle *cube = new Particle[N];
    InitPositions(cube);
    calcCM(cube,rCMStart,vCM);
    calcCM(cube,rCM,vCM);
    InitVelocities(cube);
    ComputeAccelerations(cube);
    for(run=0;run<10000;run++)
    {
            if(run%500==0)
                printf("(INIT) Zeitschritt %d\n",run);
            fprintf(positions,"%d\nTrap\n",N);
            VelocityVerlet(cube,1,positions);
            if(run%100 == 0)
                rescaleVelocities(cube);
    }
    
    for(run=0;run<10000;run++)
    {
            if(run%500==0)
                printf("(TRAP) Zeitschritt %d\n",run);
            fprintf(positions,"%d\nTrap\n",N);
            VelocityVerlet(cube,1,positions);
            //calcCM(cube,rCM,comData);
            harmonicTrap(rCM,vCM,pos,cube);
            /*
             *if(run%100 == 0)
             *{
             *    for(int n=0;n<N;n++)
             *        fprintf(positions,"Ar\t%lf\t%lf\t%lf\n",cube[n].r[0],cube[n].r[1],cube[n].r[2]);
             *}
             */
    }

    fclose(positions);


//void harmonicTrap(double* rCM, double* vCM, double* pos, Particle* particles)
	return 0;
}
