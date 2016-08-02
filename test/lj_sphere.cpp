#include <iostream>
#include <cstdio>
#include <cmath>
#include "../include/globals.hpp"
#include "../include/functions.hpp"
#include "../include/classes.hpp"

void InitPositionsNew(std::vector<Particle*> particle)
{
    unsigned int M=1;
    unsigned int n,x,y,z,k;
    n=0;
    
    while(4*M*M*M < N)
        ++M;

    double a = L/M;
    double xCell[4] = {0.25,0.75,0.75,0.25};
    double yCell[4] = {0.25,0.75,0.25,0.75};
    double zCell[4] = {0.25,0.25,0.75,0.75};
    
    for(x=0;x<M;x++)
        for(y=0;y<M;y++)
            for(z=0;z<M;z++)
                for(k=0;k<4;k++)
                    if(n<N)
                    {
                        particle[n]->r[0] = (x + xCell[k])*a;
                        particle[n]->r[1] = (y + yCell[k])*a;
                        particle[n]->r[2] = (z + zCell[k])*a;
                        particle[n]->type = 1;
                        particle[n]->name = "Glass";
						if((x==0 && (k == 3||k==0))|| (x==M-1 && (k==1 || k==2))|| (y==0 && (k==0 || k==2)) || (y==M-1 && (k==1 || k==3)) || (z==0 && (k==0 || k==1)) || (z==M-1 && (k==2 || k==3)))
							particle[n]->surface=true;
                        n++;
                    }
}

void calcCMNew(std::vector<Particle*> particles,double *rCM, double* vCM)
{
    for(unsigned int i=0;i<3;i++)
    {
        rCM[i] = 0;
    }

    for(unsigned int i=0;i<N;i++)
        for(unsigned int k=0;k<3;k++)
        {
            rCM[k]+=particles[i]->r[k];
        }

    for(unsigned int i=0;i<3;i++)
    {
        rCM[i]/=N;
    }
}

void ComputeAccelerationsNew(std::vector<Particle*> particle) 
{
    double rij[3];
    double rSqd = 0;
    for (unsigned int i = 0; i < particle.size(); i++)         // set all accelerations to zero
        for (unsigned int k = 0; k < 3; k++)
            particle[i]->a[k] = 0;
    
    for(unsigned int i=0;i<particle.size()-1;i++)
    {
        for(unsigned int j=i+1;j<particle.size();j++)
        {
            rSqd = 0;
            for(unsigned int m=0;m<3;m++)
            {
                rij[m] = particle[i]->r[m]-particle[j]->r[m];
                /*
                 *if(PBC_FLAG)
                 *{               
                 *    if(abs(rij[m]) > 0.5 * L)
                 *    {
                 *        if(rij[m] > 0)
                 *            rij[m] -= L;
                 *        else
                 *            rij[m] += L;
                 *    }
                 *}
                 */
                
                //rSqd += rij[m] * rij[m];
            }
            rSqd = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
            if(rSqd <= rCutOff*rCutOff)
            {       
                //double f = 8*48 * ( pow(rSqd,-7.) - 0.5*pow(rSqd,-4.) );
                double f = 48 * ( pow(rSqd,-7) - 0.5*pow(rSqd,-4) );
                for(unsigned int m=0; m<3; m++)
                {
                    particle[i]->a[m] += rij[m] * f;
                    particle[j]->a[m] -= rij[m] * f;
                }
            }
        }
    }
}

void VelocityVerletNew(std::vector<Particle*> particle)
{
    unsigned int i,m;

    for(i=0;i<particle.size();i++)
    {
        for(m=0;m<3;m++)
        {
            particle[i]->r[m] += particle[i]->v[m]*dt + 0.5 * particle[i]->a[m] * dt * dt;
            particle[i]->v[m] += 0.5 * particle[i]->a[m] * dt;
            
        }
    }
     
    ComputeAccelerationsNew(particle);
    for(i=0;i<particle.size();i++)
        for(m=0;m<3;m++)
            particle[i]->v[m] += 0.5 * particle[i]->a[m] * dt;
}
void InitVelocitiesNew(std::vector<Particle*> particle)
{
    unsigned int n,i;
    double vCM[3] = {0,0,0};
    double l[3] = {0,0,0};
    double rCM[3]={0,0,0};
    double omega[3]={0,0,0}; // for calculating angular velocity
    double dx[3]; // for calculating distance to rCM
    double invinertia[3][3];
    for(i=0;i<3;i++)
    {
        for(unsigned int j=0;j<3;j++)
        {
            if(i==j)
                invinertia[i][j]=(6./(L*L*particle.size()));
            else
                invinertia[i][j]=0;
        }
    }

    for(n=0;n<particle.size();n++)
        for(i=0;i<3;i++)
            particle[n]->v[i] = 2*(gsl_rng_uniform(r)-0.5); 
    
    
    for(n=0;n<particle.size();n++)
        for(i=0;i<3;i++)
            rCM[i]+=particle[n]->r[i];
    
    for(i=0;i<3;i++)
        rCM[i] /= particle.size();
        

    for(n=0;n<particle.size();n++)
        for(i=0;i<3;i++)
            vCM[i] += particle[n]->v[i];


    for(i=0;i<3;i++)
        vCM[i] /= particle.size();


    for(n=0;n<particle.size();n++)
        for(i=0;i<3;i++)
            particle[n]->v[i] -= vCM[i];

    for(n=0;n<particle.size();n++)
    {
        for(i=0;i<3;i++)
            dx[i]=rCM[i]-particle[n]->r[i];
        l[0]+= (dx[1]*particle[n]->v[2]-dx[2]*particle[n]->v[1]);
        l[1]+= (dx[2]*particle[n]->v[0]-dx[0]*particle[n]->v[2]);
        l[2]+= (dx[0]*particle[n]->v[1]-dx[1]*particle[n]->v[0]);
    }

    for(i=0;i<3;i++)
        for(unsigned int j=0;j<3;j++)
            omega[i]+=invinertia[i][j]*l[j];

    for(n=0;n<particle.size();n++)
    {
        for(i=0;i<3;i++)
            dx[i] = rCM[i] - particle[n]->r[i];
     particle[n]->v[0] -= (omega[1]*dx[2]-omega[2]*dx[1]);
     particle[n]->v[1] -= (omega[2]*dx[0]-omega[0]*dx[2]);
     particle[n]->v[2] -= (omega[0]*dx[1]-omega[1]*dx[0]);
    }

}
int main()
{
    std::vector<Particle*> particles;
    for(int i=0;i<N;i++)
        particles.push_back(new Particle());
    InitPositionsNew(particles);    
    InitVelocitiesNew(particles);
    double* center = new double[3];
    calcCMNew(particles,center,NULL);
    std::cout << "COM: (" << center[0] << ", " << center[1] << ", " << center[2] << ")" << std::endl;
    std::cout << 0.5*L << std::endl;
    double radius = 0.5*L;
    double* dist = new double[3];
    double distSqd = 0;
    for(unsigned int i=0;i<particles.size();i++)     
    {
       for(int j=0;j<3;j++)
           dist[j] = particles[i]->r[j]-center[j];
       distSqd = sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
       if(distSqd > radius)
       {
           //std::cout << i << ", ";
           if(i == 0)
               particles.erase(particles.begin());
           particles.erase(particles.begin()+i);
           i = 0;
       }
    }
    particles.erase(particles.begin());
    std::cout << std::endl;

    FILE* output = fopen("sphere_test.xyz","w");
    std::cout << "works" << std::endl; 
    ComputeAccelerationsNew(particles);
    for(int k=0;k<10000;k++)
    {
        VelocityVerletNew(particles);
        if(k%100 == 0)
        {
            fprintf(output,"%lu\nSphere\n",particles.size());
            for(unsigned int i=0;i<particles.size();i++)
                fprintf(output,"Ar\t%lf\t%lf\t%lf\n",\
                        particles[i]->r[0],particles[i]->r[1],particles[i]->r[2]);
        }
    }

    /*
     *fprintf(output,"%lu\nSphere\n",particles.size());
     *for(unsigned int i=0;i<particles.size();i++)
     *    fprintf(output,"Ar\t%lf\t%lf\t%lf\n",\
     *            particles[i]->r[0],particles[i]->r[1],particles[i]->r[2]);
     */
    fclose(output);


    return 0;
}
