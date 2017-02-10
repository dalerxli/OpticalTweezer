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

    setValues(0.2,0.04,0.2,0.8,0.06);
    Particle *particles = new Particle[N];
    InitPositions(particles);
    std::list<Particle*> gas;
    FILE* tempout = fopen("baro_temp_004_new_1_other.dat","w");
    FILE* enerout = fopen("baro_energy_004_new_1_other.dat","w");
    FILE* compos = fopen("baro_compos004_new_other.dat","w");
    FILE* comvel= fopen("baro_comvel004_new_other.dat","w");


    InitVelocities(particles);
    double* vCM = new double[3];
    double* rCM = new double[3];
    double* rCMStart = new double[3];
    double energy = 0.;
    calcCM(particles,rCMStart,vCM);
    for(int i=0;i<3;i++)
        center[i] = L/2.;

    /*
     *for(int i=0;i<200000;i++)
     *{
     *    if(i%100 == 0)
     *    {
     *        std::cout << "Step " << i << std::endl;
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
     *    if(i%100 == 0)
     *        if(i < 10000)
     *            rescaleVelocities(particles);
     *    VelocityVerlet(particles);
     *    //calcTemp(particles,tempout);
     *    //energy = calculateEnergies(particles);
     *    //fprintf(enerout,"%lf\n",energy);
     *}
     */
    

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
        //if(i%100 == 0)
            //if(i < 10000)
                //rescaleVelocities(particles);
        VelocityVerlet(particles);
        calcTemp(particles,tempout);
        energy = calculateEnergies(particles);
        fprintf(enerout,"%lf\n",energy);
    }
    
    
    verletBaroAccelerations(particles,gas);
    /*
     *for(int i=0;i<60000;i++)
     *{
     *    if(i%100 == 0)
     *    {
     *        std::cout << "Step " << i << std::endl;
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
     *    //eHEX(particles);
     *    VelocityVerlet(particles);
     *    BarostatNew(particles,gas);
     *    harmonicTrap(rCM,vCM,rCMStart,particles);
     *    calcTemp(particles,tempout);
     *    energy = calculateEnergies(particles);
     *    fprintf(enerout,"%lf\n",energy);
     *}
     */

    for(int i=0;i<60000;i++)
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
        //eHEX(particles);
        //VelocityVerlet(particles);
        //BarostatNew(particles,gas);
        //harmonicTrap(rCM,vCM,rCMStart,particles);
        eHEXBaroNewTemp(particles,gas);
        calcTemp(particles,tempout);
        energy = calculateEnergies(particles);
        fprintf(enerout,"%lf\n",energy);
    }


    fclose(tempout);
    fclose(enerout);
    fclose(compos);
    fclose(comvel);
}
