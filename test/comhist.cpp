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

    gsl_histogram2d *positionsxy = gsl_histogram2d_alloc(100,100);
    gsl_histogram2d *positionsxz = gsl_histogram2d_alloc(100,100);
    gsl_histogram2d *positionsyz = gsl_histogram2d_alloc(100,100);
    gsl_histogram2d_set_ranges_uniform(positionsxy,-2.0*L,2.0*L,-2.0*L,2.0*L);
    gsl_histogram2d_set_ranges_uniform(positionsxz,-2.0*L,2.0*L,-2.0*L,2.0*L);
    gsl_histogram2d_set_ranges_uniform(positionsyz,-2.0*L,2.0*L,-2.0*L,2.0*L);
    setValues(0.2,0.04,0.2,2,0.06);
    Particle *particles = new Particle[N];
    InitPositions(particles);
    std::list<Particle*> gas;
    FILE* tempout = fopen("comhist_temp.dat","w");
    FILE* enerout = fopen("comhist_energy.dat","w");
    FILE* compos = fopen("comhist_compos.dat","w");
    FILE* comvel= fopen("comhist_comvel.dat","w");

    InitVelocities(particles);
    double* vCM = new double[3];
    double* rCM = new double[3];
    double* vCM1 = new double[3];
    double* rCM1 = new double[3];
    double* rCMStart = new double[3];
    double energy = 0.;
    double pos[3] = {2,0,0};
    calcCM(particles,rCMStart,vCM);
    calcCM(particles,rCM1,vCM1);

    for(int i=0;i<40000;i++)
    {
        if(i%100 == 0)
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
        gsl_histogram2d_increment(positionsxy,rCM[0],rCM[1]);
        gsl_histogram2d_increment(positionsxz,rCM[0],rCM[2]);
        gsl_histogram2d_increment(positionsyz,rCM[1],rCM[2]);


        if(i%100 == 0)
            if(i < 10000)
                rescaleVelocities(particles);

        harmonicTrap(rCM1,vCM1,pos,particles);
        VelocityVerlet(particles);
        if(i >= 10000)
        {
            calcTemp(particles,tempout);
            energy = calculateEnergies(particles);
            fprintf(enerout,"%lf\n",energy);
        }
    }
    
    
    /*
     *for(int m=0;m<50;m++)
     *    InitBarostatFull(gas);
     *for(int i=0;i<30000;i++)
     *{
     *    if(i%100 == 0)
     *    {
     *        std::cout << "Step " << i << " (" << gas.size() << ")"<< std::endl;
     *        gasStatus(gas);
     *        for(int k=0;k<3;k++)
     *        {
     *            rCM[k] = 0;
     *            vCM[k] = 0;
     *        }
     *        for(unsigned int j=0;j<N;j++)
     *            for(int k=0;k<3;k++)
     *            {
     *                rCM[k] += particles[j].r[k];
     *                vCM[k] += particles[j].v[k];
     *            }
     *        for(int k=0;k<3;k++)
     *        {
     *            rCM[k] = rCM[k]/(N*1.0);
     *            vCM[k] = vCM[k]/(N*1.0);
     *        }
     *        fprintf(compos,"%lf\t%lf\t%lf\n",rCM[0],rCM[1],rCM[2]);
     *        fprintf(comvel,"%lf\t%lf\t%lf\n",vCM[0],vCM[1],vCM[2]);
     *    }
     *    if(i%100==0)
     *        GenerateOutput(particles,gas,i);
     *    //eHEX(particles);
     *    VelocityVerlet(particles);
     *    BarostatNoBoundaries(particles,gas);
     *    //harmonicTrap(rCM,vCM,rCMStart,particles);
     *    calcTemp(particles,tempout);
     *    energy = calculateEnergiesTest(particles,gas);
     *    fprintf(enerout,"%lf\n",energy);
     *}
     */
    
    FILE* comhistxy = fopen("comhistxy.dat","w");
    FILE* comhistxz = fopen("comhistxz.dat","w");
    FILE* comhistyz = fopen("comhistyz.dat","w");
    gsl_histogram2d_fprintf(comhistxy,positionsxy,"%lf","%lf");
    gsl_histogram2d_fprintf(comhistxz,positionsxz,"%lf","%lf");
    gsl_histogram2d_fprintf(comhistyz,positionsyz,"%lf","%lf");

    fclose(tempout);
    fclose(enerout);
    fclose(compos);
    fclose(comvel);
	return 0;
}
