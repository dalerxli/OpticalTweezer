#include "../include/globals.hpp"
#include "../include/classes.hpp"
#include "../include/functions.hpp"

void computeHarmonic(Particle* particle1, double* pos)
{
    double rij[3];
    double rSqd = 0;
    const double k = 15.;
    for (int k = 0; k < 3; k++)
        particle1->a[k] = 0;
    
    for(int i=0;i<3;i++)
        rij[i] = pos[i] - particle1->r[i];
    rSqd = rij[0]*rij[0]+rij[0]*rij[0]+rij[0]*rij[0];

    for(int m=0;m<3;m++)
        particle1->a[m] += k*rij[m];
}

void velocityHarmonic(Particle* particle1,double* pos)
{
    int i,m;

    for(m=0;m<3;m++)
    {
        particle1->r[m] += particle1->v[m]*dt + 0.5 * particle1->a[m] * dt * dt;
        particle1->v[m] += 0.5 * particle1->a[m] * dt;
    }
     
    computeHarmonic(particle1,pos);
    for(m=0;m<3;m++)
        particle1->v[m] += 0.5 * particle1->a[m] * dt;
}

int main()
{
    Particle *particles = new Particle[N];
    InitPositions(particles);
    FILE* positions = fopen("comtestpos.xyz","w"); 

}
