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
    setValues(0.2,0.04,0.2,2,0.06);
    Particle *particles = new Particle[N];
    InitPositions(particles);
    std::list<Particle*> gas;
    FILE* tempout = fopen("barotemp_temp.dat","w");
    FILE* enerout = fopen("barotemp_energy.dat","w");
    FILE* compos = fopen("barotemp_compos.dat","w");
    FILE* comvel= fopen("barotemp_comvel.dat","w");

    InitVelocities(particles);
    double* vCM = new double[3];
    double* rCM = new double[3];
    double* rCMStart = new double[3];
    double energy = 0.;
    calcCM(particles,rCMStart,vCM);

    for(int i=0;i<20000;i++)
    {
        if(i%100 == 0)
        {
            std::cout << "Step " << i << std::endl;
            for(int k=0;k<3;k++)
            {
                rCM[k] = 0;
                vCM[k] = 0;
            }
            for(unsigned int j=0;j<N;j++)
                for(int k=0;k<3;k++)
                {
                    rCM[k] += particles[j].r[k];
                    vCM[k] += particles[j].v[k];
                }
            for(int k=0;k<3;k++)
            {
                rCM[k] = rCM[k]/(N*1.0);
                vCM[k] = vCM[k]/(N*1.0);
            }
            fprintf(compos,"%lf\t%lf\t%lf\n",rCM[0],rCM[1],rCM[2]);
            fprintf(comvel,"%lf\t%lf\t%lf\n",vCM[0],vCM[1],vCM[2]);
        }
        if(i%100 == 0)
            if(i < 10000)
                rescaleVelocities(particles);
        VelocityVerlet(particles);
        if(i >= 10000)
        {
            calcTemp(particles,tempout);
            energy = calculateEnergies(particles);
            fprintf(enerout,"%lf\n",energy);
        }
    }
    
    
    for(int m=0;m<50;m++)
        InitBarostatFull(gas);
    for(int i=0;i<30000;i++)
    {
        if(i%100 == 0)
        {
            std::cout << "Step " << i << " (" << gas.size() << ")"<< std::endl;
            gasStatus(gas);
            for(int k=0;k<3;k++)
            {
                rCM[k] = 0;
                vCM[k] = 0;
            }
            for(unsigned int j=0;j<N;j++)
                for(int k=0;k<3;k++)
                {
                    rCM[k] += particles[j].r[k];
                    vCM[k] += particles[j].v[k];
                }
            for(int k=0;k<3;k++)
            {
                rCM[k] = rCM[k]/(N*1.0);
                vCM[k] = vCM[k]/(N*1.0);
            }
            fprintf(compos,"%lf\t%lf\t%lf\n",rCM[0],rCM[1],rCM[2]);
            fprintf(comvel,"%lf\t%lf\t%lf\n",vCM[0],vCM[1],vCM[2]);
        }
        if(i%100==0)
            GenerateOutput(particles,gas,i);
        //eHEX(particles);
        VelocityVerlet(particles);
        BarostatNoBoundaries(particles,gas);
        //harmonicTrap(rCM,vCM,rCMStart,particles);
        calcTemp(particles,tempout);
        energy = calculateEnergiesTest(particles,gas);
        fprintf(enerout,"%lf\n",energy);
    }
    
    std::cout << "\a";

    fclose(tempout);
    fclose(enerout);
    fclose(compos);
    fclose(comvel);
	return 0;
}
