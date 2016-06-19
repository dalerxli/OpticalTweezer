#include "../include/globals.hpp"
#include "../include/classes.hpp"
#include "../include/functions.hpp"

int main()
{
    Particle *particles = new Particle[N];
    InitPositions(particles);
    FILE* positions = fopen("velocitytest.xyz","w"); 
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
    for(int i=0;i<1000;i++)
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
    double* vCM = new double[3];
    double* rCM = new double[3];
    calcCM(particles,rCM,vCM);
    std::cout << "rCM: (" << rCM[0] << ", " << rCM[1] << ", " << rCM[2] << ") " << std::endl;
    std::cout << "vCM: (" << vCM[0] << ", " << vCM[1] << ", " << vCM[2] << ") " << std::endl;
    for(int i=0;i<N;i++)
        for(int k=0;k<3;k++)
            particles[i].v[k] += 0.1;
    for(int i=0;i<1000;i++)
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
    calcCM(particles,rCM,vCM);
    std::cout << "rCM: (" << rCM[0] << ", " << rCM[1] << ", " << rCM[2] << ") " << std::endl;
    std::cout << "vCM: (" << vCM[0] << ", " << vCM[1] << ", " << vCM[2] << ") " << std::endl;
    return 0;
}
