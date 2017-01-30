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
    setValues(0.2,0.04,0.2,0.5,0.06);
    Particle *particles = new Particle[N];
    InitPositions(particles);
    std::list<Particle*> gas;
    FILE* tempout = fopen("barotemp_new_temp.dat","w");
    FILE* enerout = fopen("barotemp_new_energy.dat","w");
    FILE* compos = fopen("barotemp_new_compos.dat","w");
    FILE* comvel= fopen("barotemp_new_comvel.dat","w");

    gsl_histogram2d_set_ranges_uniform(positionsxy,0,L,0,L);
    gsl_histogram2d_set_ranges_uniform(positionsxz,0,L,0,L);
    gsl_histogram2d_set_ranges_uniform(positionsyz,0,L,0,L);

    InitVelocities(particles);
    double* vCM = new double[3];
    double* rCM = new double[3];
    double* rCMStart = new double[3];
    double energy = 0.;
    calcCM(particles,rCMStart,vCM);
    for(int i=0;i<3;i++)
        center[i] = L/2.;

    for(int i=0;i<50000;i++)
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
            if(i < 10000 && i > 0)
                rescaleVelocities(particles);
        VelocityVerlet(particles);
        if(i >= 10000)
        {
            calcTemp(particles,tempout);
            energy = calculateEnergies(particles);
            fprintf(enerout,"%lf\n",energy);
        }
    }
    
    
    verletBaroAccelerations(particles,gas);
    for(int i=0;i<30000;i++)
    {
        if(i%100 == 0)
        {
            std::cout << "Step " << i << " (" << gas.size() << ")"<< std::endl;
            //gasStatus(gas);
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
        //VelocityVerlet(particles);
        //BarostatNoBoundaries(particles,gas);
        //harmonicTrap(rCM,vCM,rCMStart,particles);
        //verletBaro(particles,gas);
        eHEXBaro(particles,gas);
        calcTemp(particles,tempout);
        energy = calculateEnergiesTest(particles,gas);
        fprintf(enerout,"%lf\n",energy);
        gsl_histogram2d_increment(positionsxy,rCM[0],rCM[1]);
        gsl_histogram2d_increment(positionsxz,rCM[0],rCM[2]);
        gsl_histogram2d_increment(positionsyz,rCM[1],rCM[2]);
    }
    
    FILE* comhistxy = fopen("barotemp_new_comhistxy.dat","w");
    FILE* comhistxz = fopen("barotemp_new_comhistxz.dat","w");
    FILE* comhistyz = fopen("barotemp_new_comhistyz.dat","w");
    gsl_histogram2d_fprintf(comhistxy,positionsxy,"%lf","%lf");
    gsl_histogram2d_fprintf(comhistxz,positionsxz,"%lf","%lf");
    gsl_histogram2d_fprintf(comhistyz,positionsyz,"%lf","%lf");
    fclose(comhistxy);
    fclose(comhistxz);
    fclose(comhistyz);
    writePositions(particles,gas,"allpositions.dat");

    fclose(tempout);
    fclose(enerout);
    fclose(compos);
    fclose(comvel);
	return 0;
}
