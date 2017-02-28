#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include "../../include/globals.hpp"
#include "../../include/classes.hpp"
#include "../../include/functions.hpp"

int main()
{
    dt = 0.005;
    setValues(0.2,0.04,0.2,0.8,0.06);
    Particle *particles = new Particle[N];
    InitPositions(particles);
    FILE* tempout = fopen("verlet_test_temp.dat","w");
    FILE* enerout = fopen("verlet_test_ener.dat","w");
    FILE* compos = fopen("verlet_test_compos.dat","w");
    FILE* comvel= fopen("verlet_test_comvel.dat","w");


    InitVelocities(particles);
    double* vCM = new double[3];
    double* rCM = new double[3];
    double energy = 0.;
    //calcCM(particles,rCM,vCM);

    for(int i=0;i<200000;i++)
    {
        if(i%100 == 0)
            std::cout << "Step " << i << std::endl;
        if(i%1000 == 0)
            if(i < 20000)
                rescaleVelocities(particles);
        VelocityVerlet(particles);
        calcTemp(particles,tempout);
        energy = calculateEnergies(particles);
        fprintf(enerout,"%lf\n",energy);
    }

    fclose(tempout);
    fclose(enerout);
    fclose(compos);
    fclose(comvel);
}
