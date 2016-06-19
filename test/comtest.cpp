#include "../include/globals.hpp"
#include "../include/classes.hpp"
#include "../include/functions.hpp"

int main()
{
    Particle *particles = new Particle[N];
    InitPositions(particles);
    FILE* positions = fopen("comtestpos.xyz","w"); 
    FILE* brainfart = NULL;
    fprintf(positions,"%d\nComTest\n",N);
    for(int i=0;i<N;i++)
    {
        fprintf(positions,"Ar\t");
        for(int j=0;j<3;j++)
            fprintf(positions,"%lf\t",(particles[i].r[j])*10);
        fprintf(positions,"\n");
    }
    InitVelocities(particles);
    double* vCM = new double[3];
    calcCM(particles,rCMStart,vCM);
    for(int i=0;i<100;i++)
    {
        VelocityVerlet(particles,brainfart);
        fprintf(positions,"%d\nComTest\n",N);
        for(int i=0;i<N;i++)
        {
            fprintf(positions,"Ar\t");
            for(int j=0;j<3;j++)
                fprintf(positions,"%lf\t",(particles[i].r[j])*10);
            fprintf(positions,"\n");
        }
    }
    for(int i=0;i<N;i++)
        particles[i].r[0] += 5.5;
    for(int i=0;i<100;i++)
    {
        VelocityVerlet(particles,brainfart);
        fprintf(positions,"%d\nComTest\n",N);
        for(int i=0;i<N;i++)
        {
            fprintf(positions,"Ar\t");
            for(int j=0;j<3;j++)
                fprintf(positions,"%lf\t",(particles[i].r[j])*10);
            fprintf(positions,"\n");
        }
    }

    for(int i=0;i<10000;i++)
    {
        harmonicTrap(particles,positions);
        VelocityVerlet(particles,brainfart);
        fprintf(positions,"%d\nComTest\n",N);
        for(int k=0;k<N;k++)
        {
            fprintf(positions,"Ar\t");
            for(int j=0;j<3;j++)
                fprintf(positions,"%lf\t",(particles[k].r[j])*10);
            fprintf(positions,"\n");
        }
    }
    //harmonicTrap(particles,positions);




    return 0;
}
