#include <iostream>
#include "../include/classes.hpp"
#include "../include/globals.hpp"
#include "../include/functions.hpp"

void oneParticleVerlet(Particle* particle);
void computeOneAcceleration(Particle* particle);
//bool checkIfOnLine(Particle* particle);

int main()
{
    std::ofstream outfile;
    outfile.open("oneparticle.xyz");
    Particle* p1 = new Particle();
    for(int m=0;m<3;m++)
        p1->r[m] = 0;
    p1 -> v[0] = 0.5;
    p1 -> v[1] = 0.5;
    p1 -> v[2] = 0.5;
    for(int m=0;m<3;m++)
    {
        p1->r0[m] = p1->r[m];
        p1->v0[m] = p1->v[m];
    }

    /*
     *outfile << "1\nParticleTest\n";
     *outfile << "1\t";
     *for(int m=0;m<3;m++)
     *    outfile << p1->r[m] << " ";
     *outfile << std::endl;
     */
    
    for(int i=0;i<1000;i++)
    {
        outfile << "1\nParticleTest\n";
        outfile << "1\t";
        for(int m=0;m<3;m++)
            outfile << p1->r[m] << " ";
        outfile << std::endl;

        oneParticleVerlet(p1);
        if(!checkIfOnLine(p1))
        {
            std::cout << "Linie wurde verlassen! (" << i << ")" << std::endl;
            break;
        }
    }
    p1->v[0] = 0.1;
    for(int i=1000;i<2000;i++)
    {
        outfile << "1\nParticleTest\n";
        outfile << "1\t";
        for(int m=0;m<3;m++)
            outfile << p1->r[m] << " ";
        outfile << std::endl;

        oneParticleVerlet(p1);
        if(!checkIfOnLine(p1))
        {
            std::cout << "Linie wurde verlassen! (" << i << ")" << std::endl;
            break;
        }
    }

    return 0;
}

void oneParticleVerlet(Particle* particle)
{
    int i,m;

    for(m=0;m<3;m++)
    {
        particle->r[m] += particle->v[m]*dt + 0.5 * particle->a[m] * dt * dt;
        particle->v[m] += 0.5 * particle->a[m] * dt;
        
    }
     
    //computeOneAcceleration(particle);
    //for(i=0;i<N;i++)
        //for(m=0;m<3;m++)
            //particle->v[m] += 0.5 * particle->a[m] * dt;
}

void computeOneAcceleration(Particle* particle) 
{
    double rij[3];
    double rSqd = 0;
    for (int i = 0; i < N; i++)         // set all accelerations to zero
        for (int k = 0; k < 3; k++)
            particle[i].a[k] = 0;
    
    for(int i=0;i<N-1;i++)
    {
        for(int j=i+1;j<N;j++)
        {
            rSqd = 0;
            for(int m=0;m<3;m++)
            {
                rij[m] = particle[i].r[m]-particle[j].r[m];
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
                double f = 8*48 * ( pow(rSqd,-7.) - 0.5*pow(rSqd,-4.) );
                //double f = 48 * ( pow(rSqd,-7) - 0.5*pow(rSqd,-4) );
                for(int m=0; m<3; m++)
                {
                    particle[i].a[m] += rij[m] * f;
                    particle[j].a[m] -= rij[m] * f;
                }
            }
        }
    }
}

/*
 *bool checkIfOnLine(Particle* particle)
 *{
 *    double l;
 *    if(particle->v0[0] != 0)
 *    {
 *        l = (particle->r[0] - particle->r0[0])/particle->v0[0];
 *        if((particle->r[1] == particle->r0[1]+l*particle->v0[1]) && (particle->r[2] == particle->r0[2]+l*particle->v0[2]))
 *            return true;
 *        else
 *            return false;
 *    }
 *    else if(particle->v0[1] != 0)
 *    {
 *        l = (particle->r[1] - particle->r0[1])/particle->v0[1];
 *        if((particle->r[0] == particle->r0[0]+l*particle->v0[0]) && (particle->r[2] == particle->r0[2]+l*particle->v0[2]))
 *            return true;
 *        else
 *            return false;
 *    }
 *    else if(particle->v0[2] != 0)
 *    {
 *        l = (particle->r[2] - particle->r0[2])/particle->v0[2];
 *        if((particle->r[0] == particle->r0[0]+l*particle->v0[0]) && (particle->r[1] == particle->r0[1]+l*particle->v0[1]))
 *            return true;
 *        else
 *            return false;
 *    }
 *    else
 *        return false;
 *}
 */
