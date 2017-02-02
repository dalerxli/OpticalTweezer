#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include "../../include/functions.hpp"
#include "../../include/globals.hpp"
#include "../../include/classes.hpp"

void out(Particle* cube, int Run)
{
        std::string filename = "combined/combined_";
        filename += numberToString(Run);
        filename += ".xyz";
        /*
         *if(input != NULL)
         *{
         *    filename = "output/combined_";
         *    std::string prefix(input);
         *    filename += prefix;
         *    filename += "/combined_";
         *    filename += numberToString(Run);
         *    filename += ".xyz";
         *}
         *else
         *{
         *    std::string filename = "output/combined/combined_";
         *    filename += numberToString(Run);
         *    filename += ".xyz";
         *}
         */
        
        FILE* combinedOut = fopen(filename.c_str(),"w");
        //int printed = 0;
        fprintf(combinedOut,"%ld\nFrame 0\n",N);
        unsigned int t=0;
        for(unsigned int j=0;j<N;j++)
        {
                fprintf(combinedOut,"%s\t\%lf\t\%lf\t%lf\t%d\n",(cube[j].name).c_str(),cube[j].r[0],cube[j].r[1],cube[j].r[2],cube[j].ID);
                //printed++;
        }

        fclose(combinedOut);
}
void out(Particle* cube, Particle* gas, int Run)
{
        std::string filename = "combined/combined_";
        filename += numberToString(Run);
        filename += ".xyz";
        /*
         *if(input != NULL)
         *{
         *    filename = "output/combined_";
         *    std::string prefix(input);
         *    filename += prefix;
         *    filename += "/combined_";
         *    filename += numberToString(Run);
         *    filename += ".xyz";
         *}
         *else
         *{
         *    std::string filename = "output/combined/combined_";
         *    filename += numberToString(Run);
         *    filename += ".xyz";
         *}
         */
        
        FILE* combinedOut = fopen(filename.c_str(),"w");
        //int printed = 0;
        fprintf(combinedOut,"%ld\nFrame 0\n",N+1);
        unsigned int t=0;
        for(unsigned int j=0;j<N;j++)
        {
                fprintf(combinedOut,"%s\t\%lf\t\%lf\t%lf\t%d\n",(cube[j].name).c_str(),cube[j].r[0],cube[j].r[1],cube[j].r[2],cube[j].ID);
                //printed++;
        }
        fprintf(combinedOut,"%s\t\%lf\t\%lf\t%lf\t%d\n",
                                                (gas->name).c_str(),gas->r[0],gas->r[1],gas->r[2],gas->ID);

        fclose(combinedOut);
}
void acc(Particle* cube, Particle* gas)
{
    double fVerlet = 0;
    double fBaro = 0;
    double fHarm[3] = {0,0,0};
    double rij[3] = {0,0,0};
    double rSqd = 0;
    double rCOM[3] = {0,0,0};
    double dist[3] = {0,0,0};
    double spring = 15.;

    // INITIALIZE ACCELERATIONS

    //std::list<Particle*>::iterator gasIter;
    for(unsigned int i=0;i<N;i++)
        for(unsigned int m=0;m<3;m++)
            cube[i].aOld[m] = cube[i].a[m];
    for(unsigned int i=0;i<N;i++)
        for(unsigned int m=0;m<3;m++)
            cube[i].a[m] = 0.;        
    for(unsigned int i=0;i<3;i++)
        gas->a[i] = 0;

    // LJ INTERACTION
    for(unsigned int i=0;i<N-1;i++) 
    {
        for(unsigned int j=i+1;j<N;j++) // sum over all pairs
        {
            rSqd = 0;
            for(int m=0;m<3;m++)
              rij[m] = cube[i].r[m] - cube[j].r[m];
            rSqd = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
            if(rSqd <= rCut2)
            {
                fVerlet = 48. * (pow(rSqd,-7.) - 0.5*pow(rSqd,-4.));
                for(int m=0;m<3;m++)
                {
                    cube[i].a[m] += rij[m] * fVerlet;
                    cube[j].a[m] -= rij[m] * fVerlet;
                }
            }
        }
    }

    // SOFTSPHERE INTERACTION
    for(unsigned int i=0;i<N;i++)
    {
        rSqd = 0;
        for(int m=0;m<3;m++)
            rij[m] = gas->r[m] - cube[i].r[m];
        rSqd = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
        if(rSqd <= rCut2)
        {
            fBaro = 12. * pow(rSqd,-6.);
            for(int m=0;m<3;m++)
            {
                gas->a[m] += rij[m]*fBaro;
                cube[i].a[m] -= rij[m] * fBaro;
            }
        }

    }
}

void ver(Particle* cube, Particle* gas)
{
    for(unsigned int i=0;i<N;i++)
    {
        for(unsigned int k=0;k<3;k++)
        {
            cube[i].r[k] += cube[i].v[k]*dt + 0.5 * cube[i].a[k] * dt * dt;
            cube[i].v[k] += 0.5 * cube[i].a[k] * dt;
        }
    }
    for(unsigned int k=0;k<3;k++)
    {
            gas->r[k] += gas->v[k]*dt + 0.5 * gas->a[k] * dt * dt;
            gas->v[k] += 0.5 * gas->a[k] * dt/0.1;
    }
    acc(cube,gas);
    for(unsigned int i=0;i<N;i++)
    {
        for(unsigned int k=0;k<3;k++)
        {
            cube[i].v[k] += 0.5 * cube[i].a[k] * dt;
        }
    }
    for(unsigned int k=0;k<3;k++)
    {
            gas->v[k] += 0.5 * gas->a[k] * dt/0.1;
    }
}

void ener(Particle* cube, Particle* gas, FILE* output)
{
    double eKin = 0;
    double ePot = 0;
    double eTot = 0;
    double rSqd = 0;
    double rij[3] = {0,0,0};

    for(unsigned int i=0;i<N;i++)
    {
        rSqd = 0;
        for(int m=0;m<3;m++)
            rij[m] = gas->r[m] - cube[i].r[m];
        rSqd = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
        if(rSqd <= rCut2)
        {
            ePot += 12 * pow(rSqd,-6.);
        }
    }
    for(unsigned int k=0;k<3;k++)
        eKin += 0.1*(gas->v[0]*gas->v[0]+gas->v[1]*gas->v[1]+gas->v[2]*gas->v[2])*0.5;
    eTot = eKin + ePot;
    fprintf(output,"%lf\t%lf\t%lf\n",eKin,ePot,eTot);
}
int main(int argc, char** argv)
{
    int a;
    int run=0;
    gsl_rng_set(r,98);
    gsl_histogram_set_ranges_uniform(gas_in,-0.01,0.08);
    gsl_histogram_set_ranges_uniform(gas_out,-0.01,0.08);
    gsl_histogram_set_ranges_uniform(gas_real_in,-0.01,0.08);
    gsl_histogram2d_set_ranges_uniform(positionsxy,-2.0*L,2.0*L,-2.0*L,2.0*L);
    gsl_histogram2d_set_ranges_uniform(positionsxz,-2.0*L,2.0*L,-2.0*L,2.0*L);
    gsl_histogram2d_set_ranges_uniform(positionsyz,-2.0*L,2.0*L,-2.0*L,2.0*L);

    FILE* tempout = fopen("temperature_internal.dat","w");
    FILE* COMtempout = fopen("temperature_com.dat","w");
    FILE* comData = fopen("comdata.dat","w");
    FILE* vCOMData = fopen("vCOMData.dat","w");
    FILE* gasInData = fopen("histogram_in.dat","w");
    FILE* gasOutData = fopen("histogram_out.dat","w");
    FILE* gasRealInData = fopen("histogram_real_in.dat","w");
    FILE* gasTempData = fopen("gasTempData.dat","w");
    FILE* testgasData = fopen("testgasData.dat","w");
    Particle *cube = new Particle[N];
    double sigma = sqrt(0.0000008/0.1);
    Particle *testgas = new Particle(\
            -L,L/2.,L/2.,\
            gsl_ran_rayleigh(r,sigma),\
            gsl_ran_gaussian(r,sigma),\
            gsl_ran_gaussian(r,sigma),\
            "GasIn");
    std::list<Particle*> gas;
    std::list<Particle*> gasHistory;
    double* rCM = new double[3];
    double* rCMtemp = new double[3];
    double* vCM = new double[3];
    double* rCMharm = new double[3];
    double* vCMharm = new double[3];
    double energy = 0;
    for(int i=0;i<3;i++)
        center[i] = L/2.;

    InitPositions(cube);
    calcCM(cube,rCMStart,vCM);
    calcCM(cube,rCM,vCM);
    calcCM(cube,rCMharm,vCMharm);
    InitVelocities(cube);
    ComputeAccelerations(cube);
    for(run=0;run<20000;run++)
    {
            if(run%500==0)
                printf("(INIT) Zeitschritt %d\n",run);
            VelocityVerlet(cube,0,NULL);
            if(run%100 == 0)
            {
                if(run < 2000)
                {
                    rescaleVelocities(cube);
                }
                calcTemp(cube,tempout); 
                //calcCM(cube,rCMtemp,comData);
                /*
                 *for(int k=0;k<3;k++)
                 *{
                 *    rCM[k] = 0;
                 *    vCM[k] = 0;
                 *}
                 *for(unsigned int j=0;j<N;j++)
                 *    for(int k=0;k<3;k++)
                 *    {
                 *        rCM[k] += cube[j].r[k];
                 *        vCM[k] += cube[j].v[k];
                 *    }
                 *for(int k=0;k<3;k++)
                 *{
                 *    rCM[k] = rCM[k]/(N*1.0);
                 *    vCM[k] = vCM[k]/(N*1.0);
                 *}
                 */
                //fprintf(comData,"%lf\t%lf\t%lf\n",rCM[0],rCM[1],rCM[2]);
                //fprintf(vCOMData,"%lf\t%lf\t%lf\n",vCM[0],vCM[1],vCM[2]);
            }
            //if(run%400 == 0)
            out(cube,run);
    }


    acc(cube, testgas);
    for(run=0;run<20000;run++)
    {
            if(run%500==0)
                printf("(INIT) Zeitschritt %d\n",run);
            ver(cube,testgas);
            if(run%100 == 0)
            {
                if(run < 2000)
                {
                    rescaleVelocities(cube);
                }
                calcTemp(cube,tempout); 
                //calcCM(cube,rCMtemp,comData);
                /*
                 *for(int k=0;k<3;k++)
                 *{
                 *    rCM[k] = 0;
                 *    vCM[k] = 0;
                 *}
                 *for(unsigned int j=0;j<N;j++)
                 *    for(int k=0;k<3;k++)
                 *    {
                 *        rCM[k] += cube[j].r[k];
                 *        vCM[k] += cube[j].v[k];
                 *    }
                 *for(int k=0;k<3;k++)
                 *{
                 *    rCM[k] = rCM[k]/(N*1.0);
                 *    vCM[k] = vCM[k]/(N*1.0);
                 *}
                 */
                //fprintf(comData,"%lf\t%lf\t%lf\n",rCM[0],rCM[1],rCM[2]);
                //fprintf(vCOMData,"%lf\t%lf\t%lf\n",vCM[0],vCM[1],vCM[2]);
                ener(cube,testgas,testgasData);
            }
            //if(run%400 == 0)
            out(cube,testgas,run+20000);
    }
	return 0;
}
