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
    Particle *cube = new Particle[N];
    std::list<Particle*> gas;
    FILE* tempout = fopen("output_temp.dat","w");
    double energy = 0;
    FILE* eData = fopen("output_energies.dat","w");
    InitPositions(cube);
    FILE* initpos = fopen("initial_position.xyz","w");
    fprintf(initpos,"%d\nInitial Position\n",N);
    for(int i=0;i<N;i++)  
        fprintf(initpos,"Ar\t%lf\t%lf\t%lf\n",cube[i].r[0],cube[i].r[1],cube[i].r[2]);

    InitVelocities(cube);
    ComputeAccelerations(cube);
    for(int run=0;run<200000;run++)
    {
            if(run%500==0)
                printf("(INIT) Zeitschritt %d\n",run);
            VelocityVerlet(cube,0,NULL);
            //if(run%100 == 0 && run < 50000)
            if(run%100)
            {
                if(run < 50000)
                    rescaleVelocities(cube);
                calcTemp(cube,tempout,run); 
                energy = calculateEnergies(cube,gas); 
                fprintf(eData,"%d\t%lf\n",run,energy);
            }
            if(run%1000)
                GenerateOutput(cube,gas,run);
    }
    FILE* finalpos = fopen("final_position.xyz","w");
    fprintf(finalpos,"%d\nFinal Position\n",N);
    for(int i=0;i<N;i++)  
        fprintf(finalpos,"Ar\t%lf\t%lf\t%lf\n",cube[i].r[0],cube[i].r[1],cube[i].r[2]);

    fclose(initpos);
    fclose(finalpos);
	return 0;
}
