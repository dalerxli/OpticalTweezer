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
    setValues(0.2,0.01,0.2,0.8,0.1);
    gsl_rng_set(r,98);
    Particle *particles = new Particle[N];
    InitPositions(particles);
    std::list<Particle*> gas;
    FILE* tempout = fopen("newalgorithm.dat","w");
    FILE* gasTempData = fopen("newalgorithm_gas.dat","w");
    //FILE* enerout = fopen("baro_energy_004_trackdata_1.dat","w");
    //FILE* compos = fopen("baro_compos004_trackdata.dat","w");
    //FILE* comvel= fopen("baro_comvel004_trackdata.dat","w");
    gsl_histogram2d_set_ranges_uniform(positionsxy,0,L,0,L);
    gsl_histogram2d_set_ranges_uniform(positionsxz,0,L,0,L);
    gsl_histogram2d_set_ranges_uniform(positionsyz,0,L,0,L);
    for(int i=0;i<3;i++)
        center[i] = L/2.;

    InitVelocities(particles);
    double* vCM = new double[3];
    double* rCM = new double[3];
    //double* rCMStart = new double[3];
    //double energy = 0.;
    //calcCM(particles,rCMStart,vCM);

    ComputeAccelerations(particles);
    for(int i=0;i<10000;i++)
    {
        if(i%100 == 0)
        {
            std::cout << "Step " << i << std::endl;
            /*
             *for(int k=0;k<3;k++)
             *{
             *    rCM[k] = 0;
             *    vCM[k] = 0;
             *}
             *for(unsigned int j=0;j<N;j++)
             *    for(int k=0;k<3;k++)
             *    {
             *        rCM[k] += particles[j].r[k];
             *        vCM[k] += particles[j].v[k];
             *    }
             *for(int k=0;k<3;k++)
             *{
             *    rCM[k] = rCM[k]/(N*1.0);
             *    vCM[k] = vCM[k]/(N*1.0);
             *}
             *fprintf(compos,"%lf\t%lf\t%lf\n",rCM[0],rCM[1],rCM[2]);
             *fprintf(comvel,"%lf\t%lf\t%lf\n",vCM[0],vCM[1],vCM[2]);
             */
        }
        //if(i%400 == 0)
            //GenerateOutput(particles,gas,i);
        if(i%100 == 0)
            if(i < 10000)
                rescaleVelocities(particles);
        VelocityVerlet(particles);
        //calcTemp(particles,tempout);
        //energy = calculateEnergies(particles);
        //fprintf(enerout,"%lf\n",energy);
    }

    verletBaroAccelerations(particles,gas);
    for(int i = 0;i<10000;i++)
    {
        if(i%500==0)
        {
            printf("(MEASURE) Zeitschritt %d - Number of Gas particles: %d\n",i,gas.size());
        }
        //verletBaro(particles,gas);
        eHEXBaro(particles,gas);
        if(i%400 == 0)
            GenerateOutput(particles,gas,i);
        if(i%100==0)
        {
            calcTemp(particles,tempout); 
            //calcCOMTemp(vCM,COMtempout);
            //fprintf(vCOMData,"%lf\t%lf\t%lf\n",vCM[0],vCM[1],vCM[2]);
            calculateGasTemperature(gas,gasTempData);
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
            gsl_histogram2d_increment(positionsxy,rCM[0],rCM[1]);
            gsl_histogram2d_increment(positionsxz,rCM[0],rCM[2]);
            gsl_histogram2d_increment(positionsyz,rCM[1],rCM[2]);
        }
    }
    FILE* comhistxy = fopen("new_comhistxy.dat","w");
    FILE* comhistxz = fopen("new_comhistxz.dat","w");
    FILE* comhistyz = fopen("new_comhistyz.dat","w");
    gsl_histogram2d_fprintf(comhistxy,positionsxy,"%lf","%lf");
    gsl_histogram2d_fprintf(comhistxz,positionsxz,"%lf","%lf");
    gsl_histogram2d_fprintf(comhistyz,positionsyz,"%lf","%lf");
    fclose(comhistxy);
    fclose(comhistxz);
    fclose(comhistyz);
	return 0;
}
