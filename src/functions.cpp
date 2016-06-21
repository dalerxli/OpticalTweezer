#include "../include/functions.hpp"

void InitPositions(Particle* particle)
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
                        particle[n].r[0] = (x + xCell[k])*a;
                        particle[n].r[1] = (y + yCell[k])*a;
                        particle[n].r[2] = (z + zCell[k])*a;
                        particle[n].type = 1;
                        particle[n].name = "Glass";
						if((x==0 && (k == 3||k==0))|| (x==M-1 && (k==1 || k==2))|| (y==0 && (k==0 || k==2)) || (y==M-1 && (k==1 || k==3)) || (z==0 && (k==0 || k==1)) || (z==M-1 && (k==2 || k==3)))
							particle[n].surface=true;
                        n++;
                    }
}

unsigned int NumberOfParticles()
{
        double mu;
        //double P=0.0024;
        //This is the calculated atmospheric pressure with LJ units
        double P=5;
        double k=1.0;
        double T=0.9;
        mu = 0.01*L*L*P*pow((1/(2*M_PI*k*T)),0.5);
        //std::cout << mu << std::endl;
        return gsl_ran_poisson(r,mu);
}

void ComputeAccelerations(Particle* particle) 
{
    double rij[3];
    double rSqd = 0;
    for (unsigned int i = 0; i < N; i++)         // set all accelerations to zero
        for (unsigned int k = 0; k < 3; k++)
            particle[i].a[k] = 0;
    
    for(unsigned int i=0;i<N-1;i++)
    {
        for(unsigned int j=i+1;j<N;j++)
        {
            rSqd = 0;
            for(unsigned int m=0;m<3;m++)
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
                //double f = 8*48 * ( pow(rSqd,-7.) - 0.5*pow(rSqd,-4.) );
                double f = 48 * ( pow(rSqd,-7) - 0.5*pow(rSqd,-4) );
                for(unsigned int m=0; m<3; m++)
                {
                    particle[i].a[m] += rij[m] * f;
                    particle[j].a[m] -= rij[m] * f;
                }
            }
        }
    }
}

void VelocityVerlet(Particle* particle,int WRITE,FILE* output)
{
    unsigned int i,m;

    for(i=0;i<N;i++)
    {
        for(m=0;m<3;m++)
        {
            particle[i].r[m] += particle[i].v[m]*dt + 0.5 * particle[i].a[m] * dt * dt;
            particle[i].v[m] += 0.5 * particle[i].a[m] * dt;
            
        }
        if(WRITE==1)
        {
            fprintf(output,"1\t%lf\t%lf\t%lf\n",particle[i].r[0],particle[i].r[1],particle[i].r[2]);
        }
    }
     
    ComputeAccelerations(particle);
    for(i=0;i<N;i++)
        for(m=0;m<3;m++)
            particle[i].v[m] += 0.5 * particle[i].a[m] * dt;
}

void VelocityVerlet(Particle* particle,FILE* output)
{
    unsigned int i,m;

    for(i=0;i<N;i++)
    {
        for(m=0;m<3;m++)
        {
            particle[i].r[m] += particle[i].v[m]*dt + 0.5 * particle[i].a[m] * dt * dt;
            particle[i].v[m] += 0.5 * particle[i].a[m] * dt;
            
        }
    }
     
    ComputeAccelerations(particle);
    for(i=0;i<N;i++)
        for(m=0;m<3;m++)
            particle[i].v[m] += 0.5 * particle[i].a[m] * dt;
}

void InitVelocities(Particle* particle)
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
                invinertia[i][j]=(6./(L*L*N));
            else
                invinertia[i][j]=0;
        }
    }

    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            particle[n].v[i] = 2*(gsl_rng_uniform(r)-0.5); 
    
    
    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            rCM[i]+=particle[n].r[i];
    
    for(i=0;i<3;i++)
        rCM[i] /= N;
        

    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            vCM[i] += particle[n].v[i];


    for(i=0;i<3;i++)
        vCM[i] /= N;


    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            particle[n].v[i] -= vCM[i];

    for(n=0;n<N;n++)
    {
        for(i=0;i<3;i++)
            dx[i]=rCM[i]-particle[n].r[i];
        l[0]+= (dx[1]*particle[n].v[2]-dx[2]*particle[n].v[1]);
        l[1]+= (dx[2]*particle[n].v[0]-dx[0]*particle[n].v[2]);
        l[2]+= (dx[0]*particle[n].v[1]-dx[1]*particle[n].v[0]);
    }

    for(i=0;i<3;i++)
        for(unsigned int j=0;j<3;j++)
            omega[i]+=invinertia[i][j]*l[j];

    for(n=0;n<N;n++)
    {
        for(i=0;i<3;i++)
            dx[i] = rCM[i] - particle[n].r[i];
     particle[n].v[0] -= (omega[1]*dx[2]-omega[2]*dx[1]);
     particle[n].v[1] -= (omega[2]*dx[0]-omega[0]*dx[2]);
     particle[n].v[2] -= (omega[0]*dx[1]-omega[1]*dx[0]);
    }

}

double Pressure(Particle* particle)
{
        double sum1=0;
        double sum2=0;
        double V=N/rho;
        //double rij=0;

/*
 *    for(int i=0;i<N-1;i++)
 *    {
 *      for(int j=i+1;j<N;j++)
 *      {
 *        if(Distance(particle,i,j)<2.5*2.5)
 *        {
 *          for(int k=0;k<3;k++)
 *          {
 *            sum1+=(particle[i].r[k]-particle[j].r[k])*Force(particle,i,j)[k];     
 *          }
 *        }
 *      }
 *    }
 *
 *    return ((0.9*rho)+(1/(3*V))*sum1);
 */

    for(unsigned int i=0;i<N;i++)
        for(unsigned int j=0;j<3;j++)
        {
            sum1+=particle[i].v[j]*particle[i].v[j];
            sum2+=particle[i].r[j]*particle[i].a[j];
        }

        return ((sum1+sum2)/(3*V));     
}

double* Force(Particle *particle,int i,int j)
{
        double *f = new double[3];
        double force = 48 * ( pow(Distance(particle,i,j),-7) - 0.5*pow(Distance(particle,i,j),-4) ); 
        for(unsigned int k=0;k<3;k++)
            f[k] = particle[i].r[k]*force;

        return f;
}

double Distance(Particle* particle, int i, int j)
{
        double *r = new double[3];
        double rSqd = 0;
        for(unsigned int m=0;m<3;m++)
        {
            r[m] = particle[i].r[m]-particle[j].r[m];
            rSqd += r[m] * r[m];
        }
        return rSqd;
}

void InitBarostat(std::list<Particle*>& particles)
{
        FILE *baroOut = fopen("output/barostat.dat","a");
        //int i,j,k;
        unsigned int i;
        unsigned int N;
        double sigma=1.1;

        std::list<Particle*>::iterator iter = particles.end();

        //side: x0
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(-L/2,
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn"));
                iter++;
                //printf("%d: (%lf,%lf,%lf)\n",(*iter)->ID,(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
                fprintf(baroOut,"%lf\t%lf\t%lf\n",(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
        }

        //side xL
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(L+L/2,
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma)*(-1),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn"));
                iter++;
                //printf("%d: (%lf,%lf,%lf)\n",(*iter)->ID,(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
                fprintf(baroOut,"%lf\t%lf\t%lf\n",(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
        }   

        //side y0
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(gsl_rng_uniform(r)*L,
                            -L/2,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn"));
                iter++;
                //printf("%d: (%lf,%lf,%lf)\n",(*iter)->ID,(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
                fprintf(baroOut,"%lf\t%lf\t%lf\n",(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
        }   

        //side yL
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(gsl_rng_uniform(r)*L,
                            L+L/2,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma)*(-1),
                            gsl_ran_gaussian(r,sigma),"GasIn"));
                iter++;
                //printf("%d: (%lf,%lf,%lf)\n",(*iter)->ID,(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
                fprintf(baroOut,"%lf\t%lf\t%lf\n",(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
        }   
        
        //side z0
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,-L/2,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma),"GasIn"));
                iter++;
                //printf("%d: (%lf,%lf,%lf)\n",(*iter)->ID,(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
                fprintf(baroOut,"%lf\t%lf\t%lf\n",(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
        }   

        //side zL
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,L+L/2,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma)*(-1),"GasIn"));
                iter++;
                //printf("%d: (%lf,%lf,%lf)\n",(*iter)->ID,(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
                fprintf(baroOut,"%lf\t%lf\t%lf\n",(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
        }   
        fclose(baroOut);
        iter = particles.begin();
        for(unsigned int t=0;t<particles.size();t++)
        {
                //printf("2\t\%lf\t\%lf\t%lf\n",
                                                //(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
                iter++;
        }
}

void InitBarostat(std::vector<Particle*>& particles)
{
        //unsigned int i,j,k;
        unsigned int i;
        unsigned int N;
        double sigma=1.1;

        //side: x0
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(-L/2,
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn"));
        }

        //side xL
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(L+L/2,
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma)*(-1),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn"));
        }   

        //side y0
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(gsl_rng_uniform(r)*L,
                            -L/2,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn"));
        }   

        //side yL
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(gsl_rng_uniform(r)*L,
                            L+L/2,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma)*(-1),
                            gsl_ran_gaussian(r,sigma),"GasIn"));
        }   
        
        //side z0
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,-L/2,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma),"GasIn"));
        }   

        //side zL
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,L+L/2,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma)*(-1),"GasIn"));
        }   
}

void GenerateOutput(Particle *cube, std::list<Particle*> gas,FILE *combinedOut)
{
        // Test for File output of both cube and gas
        // goal: find a way to store all particles
        // idea:
        // -    "allocate" 2000 particles
        // -    first write out cube particles
        // -  then gas particles
        // -  give the rest some dummy position  
        const unsigned int maxParticles = 80000; 
        combinedOut = fopen("output/combined.xyz","a");
        unsigned int printed = 0;
        fprintf(combinedOut,"%d\nFrame 0\n",maxParticles);
        std::list<Particle*>::iterator gasIter = gas.begin();
        unsigned int t=0;
        for(unsigned int j=0;j<N;j++)
        {
                fprintf(combinedOut,"%s\t\%lf\t\%lf\t%lf\n",(cube[j].name).c_str(),cube[j].r[0],cube[j].r[1],cube[j].r[2]);
                printed++;
        }
        for(t=0;t<gas.size();t++)
        {
                if(printed >= maxParticles)
                {
                    std::cout << "Maximum Number of particles reached!" << std::endl;
                    break;
                }
                //printf("2\t\%lf\t\%lf\t%lf\n",
                                                //(*gasIter)->r[0],(*gasIter)->r[1],(*gasIter)->r[2]);
                fprintf(combinedOut,"%s\t\%lf\t\%lf\t%lf\n",
                                                ((*gasIter)->name).c_str(),(*gasIter)->r[0],(*gasIter)->r[1],(*gasIter)->r[2]);
                printed++;
                gasIter++;
        }
        for(unsigned int k=printed;k<maxParticles;k++)
        {
                fprintf(combinedOut,"Dummy\t%lf\t%lf\t%lf\n",L/2,L/2,L/2);
        }
        
        fclose(combinedOut);
}

void GenerateOutput(Particle *cube, std::list<Particle*> gas,int Run)
{
        //OVITO FILES
        //generate multiple files
        //with varying particle number
        //const int maxParticles = 80000; 

        std::string filename = "output/combined/combined_";
        filename += numberToString(Run);
        filename += ".xyz";
        
        FILE* combinedOut = fopen(filename.c_str(),"w");
        int printed = 0;
        fprintf(combinedOut,"%ld\nFrame 0\n",N+gas.size());
        std::list<Particle*>::iterator gasIter = gas.begin();
        unsigned int t=0;
        for(unsigned int j=0;j<N;j++)
        {
                fprintf(combinedOut,"%s\t\%lf\t\%lf\t%lf\t%d\n",(cube[j].name).c_str(),cube[j].r[0],cube[j].r[1],cube[j].r[2],cube[j].ID);
                printed++;
        }
        for(t=0;t<gas.size();t++)
        {
                //if(printed >= maxParticles)
                //{
                    //std::cout << "Maximum Number of particles reached!" << std::endl;
                    //break;
                //}
                //printf("2\t\%lf\t\%lf\t%lf\n",
                                                //(*gasIter)->r[0],(*gasIter)->r[1],(*gasIter)->r[2]);
                fprintf(combinedOut,"%s\t\%lf\t\%lf\t%lf\t%d\n",
                                                ((*gasIter)->name).c_str(),(*gasIter)->r[0],(*gasIter)->r[1],(*gasIter)->r[2],(*gasIter)->ID);

                //printed++;
                gasIter++;
        }
        //for(int k=printed;k<maxParticles;k++)
        //{
                //fprintf(combinedOut,"Dummy\t%lf\t%lf\t%lf\n",L/2,L/2,L/2);
        //}
        printDistances(cube,gas,Run);
        
        fclose(combinedOut);
}

void GenerateOutput(Particle *cube, std::vector<Particle*> gas,FILE *combinedOut)
{
        // Test for File output of both cube and gas
        // goal: find a way to store all particles
        // idea:
        // -    "allocate" 2000 particles
        // -    first write out cube particles
        // -  then gas particles
        // -  give the rest some dummy position  
        const int maxParticles = 20000; 
        combinedOut = fopen("output/combined.xyz","a");
        int printed = 0;
        fprintf(combinedOut,"%d\nFrame 0\n",maxParticles);
        unsigned int t=0;
        for(unsigned int j=0;j<N;j++)
        {
                fprintf(combinedOut,"%d\t\%lf\t\%lf\t%lf\n",cube[j].type,cube[j].r[0],cube[j].r[1],cube[j].r[2]);
                printed++;
        }
        for(t=0;t<gas.size();t++)
        {
                if(printed >= maxParticles)
                {
                    std::cout << "Maximum Number of particles reached!" << std::endl;
                    break;
                }
                //printf("2\t\%lf\t\%lf\t%lf\n",
                                                //(*gasIter)->r[0],(*gasIter)->r[1],(*gasIter)->r[2]);
                //fprintf(combinedOut,"%d\t\%lf\t\%lf\t%lf\n",
                                                //(*gasIter)->type,(*gasIter)->r[0],(*gasIter)->r[1],(*gasIter)->r[2]);

                fprintf(combinedOut,"%d\t\%lf\t\%lf\t%lf\n",
                                                gas[t]->type,gas[t]->r[0],gas[t]->r[1],gas[t]->r[2]);
                printed++;
        }
        for(unsigned int k=printed;k<maxParticles;k++)
                fprintf(combinedOut,"8\t%lf\t%lf\t%lf\n",L/2,L/2,L/2);
        
        fclose(combinedOut);
}

void CheckBoundaries(std::list<Particle*> &particles)
{
    std::list<Particle*>::iterator temp;
    double lowerBound = -(L/2 + eps);
    double upperBound = L + L/2 + eps;

    for(temp=particles.begin();temp!=particles.end();)
    {
        if((*temp)->r[0] < lowerBound || (*temp)->r[0] > upperBound || 
                (*temp)->r[1] < lowerBound || (*temp)->r[1] > upperBound 
                ||(*temp)->r[2] < lowerBound || (*temp)->r[2] > upperBound)
        {
            //std::cout << "Deleted [ID]: " << (*temp)->ID << std::endl;
            temp=particles.erase(temp);
        }
        else
            ++temp;
    }

}

void CheckBoundaries(Particle* cube,std::list<Particle*> &particles)
{
    std::list<Particle*>::iterator temp;
    double lowerBound = -(L/2 + eps);
    double upperBound = L + L/2 + eps;
    //double dist[3] = {0,0,0};
    //double distSqd = 0.;
    //double vSqd = 0.;

    for(temp=particles.begin();temp!=particles.end();)
    {
        //if(!checkIfOnLine(temp))
            //(*temp)->name = "GasOut";
        /*
         *vSqd = 0.;
         *for(int m=0;m<3;m++)
         *    vSqd += (*temp)->v[m];
         *if(vSqd < 1e-4)
         *    (*temp)->name = "GasOut";
         */


        /*
         *for(int i=0;i<N;i++) // iterate over all particles
         *{
         *    if(cube[i].surface && std::strcmp(((*temp)->name).c_str(),"GasOut") != 0) // if particle is on the surface
         *    {
         *       for(int m=0;m<3;m++)
         *          dist[m] = (*temp)->r[m] - cube[i].r[m];
         *      distSqd = dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2];
         *      if(distSqd < 1.2)
         *      {
         *          (*temp)->name = "GasOut";
         *          //std::cout << (*temp)->ID << std::endl;
         *      }
         *    }
         *}
         */
        if((*temp)->r[0] < lowerBound || (*temp)->r[0] > upperBound || 
                (*temp)->r[1] < lowerBound || (*temp)->r[1] > upperBound 
                ||(*temp)->r[2] < lowerBound || (*temp)->r[2] > upperBound)
        {
            //std::cout << "Deleted [ID]: " << (*temp)->ID << std::endl;
            temp=particles.erase(temp);
        }
        else
            ++temp;
    }

}

void printDistances(Particle* cube, std::list<Particle*> &particles,int Run)
{
    std::list<Particle*>::iterator temp;
    double dist[3] = {0,0,0};
    double distSqd = 0.;
    std::string filename = "output/combined/distances_";
    filename += numberToString(Run);
    filename += ".dat";
    FILE* output = fopen(filename.c_str(),"w");
    
    for(temp=particles.begin();temp!=particles.end();temp++)
    {
        for(unsigned int i=0;i<N;i++)
        {
            if(cube[i].surface)
            {
                for(unsigned int m=0;m<3;m++)
                    dist[m] = (*temp)->r[m] - cube[i].r[m];
                distSqd = dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2];
                if(distSqd < 2)
                {
                    //fprintf(output,"Gas %d\tCube %d\t%lf\n",(*temp)->ID,i,distSqd); 
                    fprintf(output,"%lf\n",distSqd); 
                }
            }
        }
    }
    fclose(output);


}

void CheckBoundaries(std::vector<Particle*> &particles)
{
    double lowerBound = -(L/2 + eps);
    double upperBound = L + L/2 + eps;

    /*
     *for(temp=particles.begin();temp!=particles.end();)
     *{
     *    if((*temp)->r[0] < lowerBound || (*temp)->r[0] > upperBound || 
     *            (*temp)->r[1] < lowerBound || (*temp)->r[1] > upperBound 
     *            ||(*temp)->r[2] < lowerBound || (*temp)->r[2] > upperBound)
     *    {
     *        //std::cout << "Deleted [ID]: " << (*temp)->ID << std::endl;
     *        temp=particles.erase(temp);
     *    }
     *    else
     *        ++temp;
     *}
     */

    for(unsigned int i=0;i<particles.size();i++)
    {
        if(particles[i]->r[0] < lowerBound || particles[i]->r[0] > upperBound || 
                 particles[i]->r[1] < lowerBound || particles[i]->r[1] > upperBound 
                 ||particles[i]->r[2] < lowerBound || particles[i]->r[2] > upperBound)
        {
            particles.erase(particles.begin()+i);
        }

    }
        
}

void Barostat(Particle* cube, std::list<Particle*>& gas)
{
    std::list<Particle*>::iterator gasIter;
    InitBarostat(gas);
    
    for(gasIter = gas.begin();gasIter != gas.end(); gasIter++)
    {
        for(unsigned int m=0;m<3;m++)
        {
            (*gasIter)->r[m] += (*gasIter)->v[m]*dt + 0.5 * (*gasIter)->a[m] * dt * dt;
            (*gasIter)->v[m] += 0.5 * (*gasIter)->a[m] * dt;
        }
    }

    ComputeSoftSphere(gas,cube);
    for(gasIter = gas.begin();gasIter != gas.end(); gasIter++)
        for(unsigned int m=0;m<3;m++)
            (*gasIter)->v[m] += 0.5 * (*gasIter)->a[m] * dt;
    CheckBoundaries(cube,gas);
}

void ComputeSoftSphere(std::list<Particle*>& gas, Particle* cube)
{
    std::list<Particle*>::iterator iter;
    double rij[3];
    double rSqd=0;
    double f;
    for (iter=gas.begin();iter!=gas.end();iter++)         // set all accelerations to zero
    {
        for(unsigned int i=0;i<3;i++)
        {
            (*iter)->a[i]=0;
            (*iter)->type = 3;
        }
    }

    
    for (iter=gas.begin();iter!=gas.end();iter++)
    {
        (*iter)->potE = 0.;
        for(unsigned int n=0;n<N;n++)
        {
            rSqd=0;
            for(unsigned int i=0;i<3;i++)
            {
                rij[i] = (*iter)->r[i] - cube[n].r[i]; 
                rSqd += rij[i]*rij[i];
            }
            if(rSqd < 2.5*2.5)
            {
                (*iter)->potE += pow(rSqd,-6);
                (*iter)->type = 4;
                f = 12 * pow(rSqd,-6);
                for(unsigned int m=0;m<3;m++)
                {
                    (*iter) -> a[m] += rij[m]*f;
                    cube[n].a[m] -= rij[m] *f;
                }
            }
        }
    }
}

void PrintAllData(Particle* cube, std::list<Particle*> gas,FILE* output)
{
    for(unsigned int i=0;i<N;i++)
    {
        fprintf(output,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",i,
                cube[i].r[0],cube[i].r[1],cube[i].r[2],
                cube[i].v[0],cube[i].v[1],cube[i].v[2],
                cube[i].a[0],cube[i].a[1],cube[i].a[2]);
    }

    for(std::list<Particle*>::iterator iter=gas.begin();iter != gas.end();iter++)
    {
        fprintf(output,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",(*iter)->ID,
                (*iter)->r[0],(*iter)->r[1],(*iter)->r[2],
                (*iter)->v[0],(*iter)->v[1],(*iter)->v[2],
                (*iter)->a[0],(*iter)->a[1],(*iter)->a[2]);
    }

}

void eHEX(Particle* cube,FILE* output)
{
	double eps[N][3];
	double K=0;
	double eta[N][3];
	double xi=0;
	double sumforces[3]={0,0,0};
	double sumvel=0; 
	unsigned int n,i; 


	for(n=0;n<N;n++)
		for(i=0;i<3;i++)
			K+=(cube[n].v[i]*cube[n].v[i]);
	K=K/2;
	if(K==0)
		printf("Error, K must not be 0!");
	xi = sqrt(1+dQ/K);

	for(n=0;n<N;n++)
		for(i=0;i<3;i++)
			eta[n][i]=(FG/(2*K))*cube[n].v[i];

	for(n=0;n<N;n++)
		for(i=0;i<3;i++)
		{
			sumforces[i]+=cube[n].a[i];
			sumvel+=cube[n].v[i]*cube[n].a[i];
		}

	for(n=0;n<N;n++)
	{
		for(i=0;i<3;i++)
		{
			eps[n][i]=(FG/(12*K))*(1-(cube[n].v[i]*cube[n].v[i]/K))*
                (cube[n].a[i]-sumforces[i]/N+eta[n][i]/4*(1+(3./N)));
			//eps[n][i]=(eta[n][i]/K)*((FG/48)+(1/6)*sumvel)-(FG/(12*K))*
            //(cube[n].a[i]-(sumforces[i]));
		}
	}

	//on first run, calculate rnew_0 and vnew_0
	if(!EHEX_FLAG)
	{
		for(n=0;n<N;n++)
		{
			for(i=0;i<3;i++)
			{
				cube[n].rnew[i]=cube[n].r[i]-dt*dt*dt*eps[n][i];
				cube[n].vnew[i]=xi*cube[n].v[i]; //v_gamma=0 
			}
		}
		//printf("\nFirst run of eHEX\n");
		EHEX_FLAG=true;
	}	

	for(n=0;n<N;n++)
	{
		for(i=0;i<3;i++)
		{
			cube[n].vhalf[i]=cube[n].vnew[i]+dt*0.5*cube[n].a[i];
			cube[n].r[i]=cube[n].rnew[i]+dt*cube[n].vhalf[i];
		}	
        //fprintf(output,"cube %lf %lf %lf\n",cube[n].r[0],cube[n].r[1],cube[n].r[2]);
	}
	
    ComputeAccelerations(cube);

	for(n=0;n<N;n++)
		for(i=0;i<3;i++)
			cube[n].v[i] = cube[n].vhalf[i]+dt*0.5*cube[n].a[i];
	
	for(n=0;n<N;n++)
		for(i=0;i<3;i++)
		{
			cube[n].rnew[i]=cube[n].r[i]-dt*dt*dt*eps[n][i];
			cube[n].vnew[i]=xi*cube[n].v[i]; //v_gamma=0 
		}
}

void calcTemp(Particle* cube,FILE* output)
{
    double T=0;
    for(unsigned int i=0;i<N;i++)
        for(unsigned int k=0;k<3;k++)
            T+=cube[i].v[k]*cube[i].v[k];
    T = T/(3*N);
    fprintf(output,"%lf\n",T);
}

void calcCM(Particle* particles,double *rCM, double* vCM)
{
    for(unsigned int i=0;i<3;i++)
    {
        rCM[i] = 0;
        vCM[i] = 0;
    }

    for(unsigned int i=0;i<N;i++)
        for(unsigned int k=0;k<3;k++)
        {
            rCM[k]+=particles[i].r[k];
            vCM[k]+=particles[i].v[k];
        }

    for(unsigned int i=0;i<3;i++)
    {
        rCM[i]/=N;
        vCM[i]/=N;
    }
}

bool fileExist(const std::string& filename)
{
    struct stat buffer;
    return(stat(filename.c_str(),&buffer) == 0);
}

double totErg(Particle* particles)
{
    double ener1=0;
    double ener2=0;
    double r=0;

    for(unsigned int i=0;i<N;i++)
        for(unsigned int j=0;j<3;j++)
        {
            ener1+=particles[i].v[j]*particles[i].v[j];
        }
    ener1/=2;
    for(unsigned int i=0;i<N-1;i++)
    {
        for(unsigned int j=i+1;j<N;j++)
        {
            r = 0;
            for(unsigned int k=0;k<3;k++)
                r+=(particles[i].r[k]-particles[j].r[k])*(particles[i].r[k]-particles[j].r[k]);
            r = sqrt(r);
            ener2+=4*(pow(r,-12)-pow(r,-6));
        }
    }

    return ener1+ener2;
}

void writePositions(Particle* particles, std::string filename)
{
    FILE* writeFile = fopen(filename.c_str(),"wb");
    for(unsigned int i=0;i<N;i++)
    {
        for(unsigned int j=0;j<3;j++)
            fwrite(&(particles[i].r[j]),sizeof(double),1,writeFile);
        for(unsigned int j=0;j<3;j++)
            fwrite(&(particles[i].v[j]),sizeof(double),1,writeFile);
        for(unsigned int j=0;j<3;j++)
            fwrite(&(particles[i].a[j]),sizeof(double),1,writeFile);
        fwrite(&(particles[i].surface),sizeof(bool),1,writeFile);
    }
    fclose(writeFile);
}

void readPositions(Particle* particles, std::string filename)
{
    FILE* readFile = fopen(filename.c_str(),"rb");
    for(unsigned int i=0;i<N;i++)
    {
        for(unsigned int j=0;j<3;j++)
            fread(&(particles[i].r[j]),sizeof(double),1,readFile);
        for(unsigned int j=0;j<3;j++)
            fread(&(particles[i].v[j]),sizeof(double),1,readFile);
        for(unsigned int j=0;j<3;j++)
            fread(&(particles[i].a[j]),sizeof(double),1,readFile);
        fread(&(particles[i].surface),sizeof(bool),1,readFile);
    }
    fclose(readFile);
}

void harmonicTrap(double* rCM, double* vCM, double* pos, Particle* particles)
{
    double rij[3];
    //double rSqd = 0;
    double rCMtemp[3] = {0,0,0};
    double aCM[3] = {0,0,0};
    const double k = 5.;

    for(unsigned int i=0;i<3;i++)
        rCMtemp[i] = rCM[i];
    
    
    for(unsigned int i=0;i<3;i++)
        rij[i] = pos[i] - rCM[i];

    for(unsigned int m=0;m<3;m++)
        aCM[m] += k*rij[m];

    for(unsigned int m=0;m<3;m++)
    {
        rCM[m] += vCM[m]*dt + 0.5 * aCM[m] * dt * dt;
        vCM[m] += 0.5 * aCM[m] * dt;
    }
     
    for(unsigned int i=0;i<3;i++)
        rij[i] = pos[i] - rCM[i];
    //rSqd = rij[0]*rij[0]+rij[0]*rij[0]+rij[0]*rij[0];

    aCM[0] = 0;
    aCM[1] = 0;
    aCM[2] = 0;
    for(unsigned int m=0;m<3;m++)
        aCM[m] += k*rij[m];

    for(unsigned int m=0;m<3;m++)
        vCM[m] += 0.5 * aCM[m] * dt;

    for(unsigned int i=0;i<N;i++)
     for(unsigned int k=0;k<3;k++)
         particles[i].r[k] += (rCMtemp[k] - rCM[k]);
}

void trackParticle(Particle* cube, std::list<Particle*> gas, int partID, FILE* output)
{
    /*
     *  VERSION 1 -- Cube and Gas particles
     *std::list<Particle*>::iterator gasIter;
     *for(gasIter = gas.begin();gasIter!=gas.end();gasIter++)
     *{
     *    if((*gasIter)->ID == partID)
     *    {
     *        fprintf(output,"%d\ntest\n",N+1);
     *        for(unsigned int j=0;j<N;j++)
     *        {
     *                fprintf(output,"%s\t\%lf\t\%lf\t%lf\n",(cube[j].name).c_str(),cube[j].r[0],cube[j].r[1],cube[j].r[2]);
     *        }
     *        fprintf(output,"%s\t\%lf\t\%lf\t%lf\n",
     *                ((*gasIter)->name).c_str(),(*gasIter)->r[0],(*gasIter)->r[1],(*gasIter)->r[2]);
     *    }
     *}
     */

     std::list<Particle*>::iterator gasIter;
     for(gasIter = gas.begin();gasIter!=gas.end();gasIter++)
     {
         if((*gasIter)->ID == partID)
         {
             fprintf(output,"1\ntest\n");
             /*
              *fprintf(output,"%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
              *        ((*gasIter)->name).c_str(),(*gasIter)->r[0],(*gasIter)->r[1],(*gasIter)->r[2],\
              *        (*gasIter)->v[0],(*gasIter)->v[1],(*gasIter)->v[2],sqrt((*gasIter)->v[0]*(*gasIter)->v[0]+(*gasIter)->v[1]*(*gasIter)->v[1]+(*gasIter)->v[2]*(*gasIter)->v[2]),\
              *        (*gasIter)->a[0],(*gasIter)->a[1],(*gasIter)->a[2],sqrt((*gasIter)->a[0]*(*gasIter)->a[0]+(*gasIter)->a[1]*(*gasIter)->a[1]+(*gasIter)->a[2]*(*gasIter)->a[2]),\
              *        (*gasIter)->potE),(*gasIter)->m*((*gasIter)->v[0]*(*gasIter)->v[0]+(*gasIter)->v[1]*(*gasIter)->v[1]+(*gasIter)->v[2]*(*gasIter)->v[2])/2.;
              */

             fprintf(output,"%s",((*gasIter)->name).c_str());
             fprintf(output,"\t%lf\t%lf\t%lf",(*gasIter)->r[0],(*gasIter)->r[1],(*gasIter)->r[2]);
             fprintf(output,"\t%lf\t%lf\t%lf",(*gasIter)->v[0],(*gasIter)->v[1],(*gasIter)->v[2]);
             fprintf(output,"\t%lf",sqrt((*gasIter)->v[0]*(*gasIter)->v[0]+(*gasIter)->v[1]*(*gasIter)->v[1]+(*gasIter)->v[2]*(*gasIter)->v[2]));
             fprintf(output,"\t%lf\t%lf\t%lf",(*gasIter)->a[0],(*gasIter)->a[1],(*gasIter)->a[2]);
             fprintf(output,"\t%lf",sqrt((*gasIter)->a[0]*(*gasIter)->a[0]+(*gasIter)->a[1]*(*gasIter)->a[1]+(*gasIter)->a[2]*(*gasIter)->a[2]));
             fprintf(output,"\t%lf",(*gasIter)->potE);
             fprintf(output,"\t%lf",(*gasIter)->m*((*gasIter)->v[0]*(*gasIter)->v[0]+(*gasIter)->v[1]*(*gasIter)->v[1]+(*gasIter)->v[2]*(*gasIter)->v[2])/2.);
             fprintf(output,"\n");


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


bool checkIfOnLine(std::list<Particle*>::iterator iterator)
{
    double l;
    if((*iterator)->v0[0] != 0)
    {
        l = ((*iterator)->r[0] - (*iterator)->r0[0])/(*iterator)->v0[0];
        if(((*iterator)->r[1] == (*iterator)->r0[1]+l*(*iterator)->v0[1]) && ((*iterator)->r[2] == (*iterator)->r0[2]+l*(*iterator)->v0[2]))
            return true;
        else
            return false;
    }
    else if((*iterator)->v0[1] != 0)
    {
        l = ((*iterator)->r[1] - (*iterator)->r0[1])/(*iterator)->v0[1];
        if(((*iterator)->r[0] == (*iterator)->r0[0]+l*(*iterator)->v0[0]) && ((*iterator)->r[2] == (*iterator)->r0[2]+l*(*iterator)->v0[2]))
            return true;
        else
            return false;
    }
    else if((*iterator)->v0[2] != 0)
    {
        l = ((*iterator)->r[2] - (*iterator)->r0[2])/(*iterator)->v0[2];
        if(((*iterator)->r[0] == (*iterator)->r0[0]+l*(*iterator)->v0[0]) && ((*iterator)->r[1] == (*iterator)->r0[1]+l*(*iterator)->v0[1]))
            return true;
        else
            return false;
    }
    else
        return false;
}
