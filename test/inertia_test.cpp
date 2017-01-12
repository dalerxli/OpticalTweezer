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

int main(int argc, char** argv)
{
    int run=0;
    Particle *cube = new Particle[N];
    InitPositions(cube);
    InitVelocities(cube);
    ComputeAccelerations(cube);
    FILE* inertiatest = fopen("intertia_test.xyz","w");
    fprintf(inertiatest,"%d\nCOMMENT\n",N);
    for(int i=0;i<N;i++)
        fprintf(inertiatest,"Ar\t%lf\t%lf\t%lf\n",cube[i].r[0],cube[i].r[1],cube[i].r[2]);

    for(run=0;run<100000;run++)
    {
            if(run%500==0)
                printf("(INIT) Zeitschritt %d\n",run);
            VelocityVerlet(cube,0,NULL);
            if(run%100 == 0)
            {
                //rescaleVelocities(cube);
                fprintf(inertiatest,"%d\nCOMMENT\n",N);
                for(int i=0;i<N;i++)
                    fprintf(inertiatest,"Ar\t%lf\t%lf\t%lf\n",cube[i].r[0],cube[i].r[1],cube[i].r[2]);
            }
    }
    fprintf(inertiatest,"%d\nCOMMENT\n",N);
    for(int i=0;i<N;i++)
        fprintf(inertiatest,"Ar\t%lf\t%lf\t%lf\n",cube[i].r[0],cube[i].r[1],cube[i].r[2]);
	return 0;
}
