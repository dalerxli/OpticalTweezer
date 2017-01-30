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
        //double P=4;
        double k=1.0;
        double m = 0.1;
        //double T=0.9;
        mu = dt*3*L*3*L*P*pow((1./(2*m*M_PI*k*AmbientTemp)),0.5);
        //std::cout << mu << std::endl;
        return gsl_ran_poisson(r,mu);
}

void ComputeAccelerations(Particle* particle) 
{
    double f = 0;
    double v = 0;
    double rij[3];
    double rSqd = 0;
    for (unsigned int i = 0; i < N; i++)         // set all accelerations to zero
        for (unsigned int k = 0; k < 3; k++)
            particle[i].a[k] = 0;
    /*
     *for(unsigned int i=0;i<N;i++)
     *    for(unsigned int j=0;j<N;j++)
     *        for(unsigned int k=0;k<3;k++)
     *            Forces[i][j]
     */
    
    for(unsigned int i=0;i<N-1;i++)
    {
        for(unsigned int j=i+1;j<N;j++)
        {
            rSqd = 0;
            for(unsigned int m=0;m<3;m++)
            {
                rij[m] = particle[i].r[m]-particle[j].r[m];
                //Distances[i][j][m] = rij[m];
                //Distances[j][i][m] = (-1.)*rij[m];
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
                //f = 8*48 * ( pow(rSqd,-7.) - 0.5*pow(rSqd,-4.) );
                f = 48 * ( pow(rSqd,-7) - 0.5*pow(rSqd,-4) );
                for(unsigned int m=0; m<3; m++)
                {
                    particle[i].a[m] += rij[m] * f;
                    particle[j].a[m] -= rij[m] * f;
                }
            }
        }
    }

    /*
     *for(unsigned int i=0;i<N-1;i++)
     *    for(unsigned int j=i+i;j<N;j++)
     *        for(unsigned int m=0;m<3;m++)
     *        {
     *                //Forces[i][j][m] = particle[i].a[m];
     *                //Forces[j][i][m] = particle[j].a[m];
     *        }
     */
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

void VelocityVerlet(Particle* particle)
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
    /*
     *for(i=0;i<3;i++)
     *{
     *    for(unsigned int j=0;j<3;j++)
     *    {
     *        if(i==j)
     *            invinertia[i][j]=(6./(L*L*N));
     *        else
     *            invinertia[i][j]=0;
     *    }
     *}
     */

    /*
     *for(n=0;n<N;n++)
     *    for(i=0;i<3;i++)
     *        particle[n].v[i] = 2*(gsl_rng_uniform(r)-0.5); 
     */
    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
           particle[n].v[i] = gsl_ran_gaussian(r,0.01)+Temp;
    
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

    /*
     *for(n=0;n<N;n++)
     *{
     *    for(i=0;i<3;i++)
     *        dx[i]=rCM[i]-particle[n].r[i];
     *    l[0]+= (dx[1]*particle[n].v[2]-dx[2]*particle[n].v[1]);
     *    l[1]+= (dx[2]*particle[n].v[0]-dx[0]*particle[n].v[2]);
     *    l[2]+= (dx[0]*particle[n].v[1]-dx[1]*particle[n].v[0]);
     *}
     */

/*
 *    for(i=0;i<3;i++)
 *        for(unsigned int j=0;j<3;j++)
 *            omega[i]+=invinertia[i][j]*l[j];
 *
 *    for(n=0;n<N;n++)
 *    {
 *        for(i=0;i<3;i++)
 *            dx[i] = rCM[i] - particle[n].r[i];
 *     particle[n].v[0] -= (omega[1]*dx[2]-omega[2]*dx[1]);
 *     particle[n].v[1] -= (omega[2]*dx[0]-omega[0]*dx[2]);
 *     particle[n].v[2] -= (omega[0]*dx[1]-omega[1]*dx[0]);
 *    }
 */

    //rescaleVelocities(particle);

}

void InitVelocitiesTest(Particle* particle)
{
    unsigned int n,i;
    double vCM[3] = {0,0,0};
    double l[3] = {0,0,0};
    double rCM[3]={0,0,0};
    double omega[3]={0,0,0}; // for calculating angular velocity
    double** om = new double*[N];
    for(n=0;n<N;n++)
        om[n] = new double[3];
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
        om[n][0]+= (dx[1]*particle[n].v[2]-dx[2]*particle[n].v[1]);
        om[n][1]+= (dx[2]*particle[n].v[0]-dx[0]*particle[n].v[2]);
        om[n][2]+= (dx[0]*particle[n].v[1]-dx[1]*particle[n].v[0]);
        for(i=0;i<3;i++)
            om[n][i] = om[n][i]/(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
    }


/*
 *    for(n=0;n<N;n++)
 *    {
 *        for(i=0;i<3;i++)
 *            dx[i]=rCM[i]-particle[n].r[i];
 *        l[0]+= (dx[1]*particle[n].v[2]-dx[2]*particle[n].v[1]);
 *        l[1]+= (dx[2]*particle[n].v[0]-dx[0]*particle[n].v[2]);
 *        l[2]+= (dx[0]*particle[n].v[1]-dx[1]*particle[n].v[0]);
 *    }
 *    //for(i=0;i<3;i++)
 *        //l[i] = l[i]/(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
 *
 *    for(i=0;i<3;i++)
 *        for(unsigned int j=0;j<3;j++)
 *            omega[i]+=invinertia[i][j]*l[j];
 *    for(i=0;i<3;i++)
 *        omega[i] = omega[i]/(N*1.0);
 */

    for(n=0;n<N;n++)
    {
        for(i=0;i<3;i++)
            dx[i] = rCM[i] - particle[n].r[i];
         particle[n].v[0] -= (om[n][1]*dx[2]-om[n][2]*dx[1]);
         particle[n].v[1] -= (om[n][2]*dx[0]-om[n][0]*dx[2]);
         particle[n].v[2] -= (om[n][0]*dx[1]-om[n][1]*dx[0]);
    }

    for(i=0;i<3;i++)
        vCM[i] = 0;

    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            vCM[i] += particle[n].v[i];


    for(i=0;i<3;i++)
        vCM[i] /= N;


    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            particle[n].v[i] -= vCM[i];

    rescaleVelocities(particle);

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
        //FILE *baroOut = fopen("output/barostat.dat","a");
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

                /*
                 *particles.push_back(new Particle(-L/2.,
                 *            (2*(gsl_rng_uniform(r)-0.5)*L)+0.5*L,
                 *            (2*(gsl_rng_uniform(r)-0.5)*L)+0.5*L,
                 *            gsl_ran_rayleigh(r,sigma),
                 *            gsl_ran_gaussian(r,sigma),
                 *            gsl_ran_gaussian(r,sigma),"GasIn"));
                 */
                iter++;
                //printf("%d: (%lf,%lf,%lf)\n",(*iter)->ID,(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
                //fprintf(baroOut,"%lf\t%lf\t%lf\n",(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
        }

        //side xL
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(L+L/2.,
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma)*(-1.),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn"));
                /*
                 *particles.push_back(new Particle(L+L/2.,
                 *            (2*(gsl_rng_uniform(r)-0.5)*L)+0.5*L,
                 *            (2*(gsl_rng_uniform(r)-0.5)*L)+0.5*L,
                 *            gsl_ran_rayleigh(r,sigma)*(-1.),
                 *            gsl_ran_gaussian(r,sigma),
                 *            gsl_ran_gaussian(r,sigma),"GasIn"));
                 */
                iter++;
                //printf("%d: (%lf,%lf,%lf)\n",(*iter)->ID,(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
                //fprintf(baroOut,"%lf\t%lf\t%lf\n",(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
        }   

        //side y0
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(gsl_rng_uniform(r)*L,
                            -L/2.,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn"));
                /*
                 *particles.push_back(new Particle((2*(gsl_rng_uniform(r)-0.5)*L)+0.5*L,
                 *            -L/2.,
                 *            (2*(gsl_rng_uniform(r)-0.5)*L)+0.5*L,
                 *            gsl_ran_gaussian(r,sigma),
                 *            gsl_ran_rayleigh(r,sigma),
                 *            gsl_ran_gaussian(r,sigma),"GasIn"));
                 */
                iter++;
                //printf("%d: (%lf,%lf,%lf)\n",(*iter)->ID,(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
                //fprintf(baroOut,"%lf\t%lf\t%lf\n",(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
        }   

        //side yL
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(gsl_rng_uniform(r)*L,
                            L+L/2.,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma)*(-1.),
                            gsl_ran_gaussian(r,sigma),"GasIn"));
                /*
                 *particles.push_back(new Particle((2*(gsl_rng_uniform(r)-0.5)*L)+0.5*L,
                 *            L+L/2.,
                 *            (2*(gsl_rng_uniform(r)-0.5)*L)+0.5*L,
                 *            gsl_ran_gaussian(r,sigma),
                 *            gsl_ran_rayleigh(r,sigma)*(-1.),
                 *            gsl_ran_gaussian(r,sigma),"GasIn"));
                 */
                iter++;
                //printf("%d: (%lf,%lf,%lf)\n",(*iter)->ID,(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
                //fprintf(baroOut,"%lf\t%lf\t%lf\n",(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
        }   
        
        //side z0
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            -L/2.,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma),"GasIn"));
                /*
                 *particles.push_back(new Particle((2*(gsl_rng_uniform(r)-0.5)*L)+0.5*L,
                 *            (2*(gsl_rng_uniform(r)-0.5)*L)+0.5*L,
                 *            -L/2.,
                 *            gsl_ran_gaussian(r,sigma),
                 *            gsl_ran_gaussian(r,sigma),
                 *            gsl_ran_rayleigh(r,sigma),"GasIn"));
                 */
                iter++;
                //printf("%d: (%lf,%lf,%lf)\n",(*iter)->ID,(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
                //fprintf(baroOut,"%lf\t%lf\t%lf\n",(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
        }   

        //side zL
        N = NumberOfParticles();
        //std::cout << "N = " << N << std::endl;
        for(i=0;i<N;i++)
        {
                particles.push_back(new Particle(
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            L+L/2.,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma)*(-1.),"GasIn"));
                /*
                 *particles.push_back(new Particle((2*(gsl_rng_uniform(r)-0.5)*L)+0.5*L,
                 *            (2*(gsl_rng_uniform(r)-0.5)*L)+0.5*L,
                 *            L+L/2.,
                 *            gsl_ran_gaussian(r,sigma),
                 *            gsl_ran_gaussian(r,sigma),
                 *            gsl_ran_rayleigh(r,sigma)*(-1.),"GasIn"));
                 */
                iter++;
                //printf("%d: (%lf,%lf,%lf)\n",(*iter)->ID,(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
                //fprintf(baroOut,"%lf\t%lf\t%lf\n",(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
        }   
        //fclose(baroOut);
        //iter = particles.begin();
        /*
         *for(unsigned int t=0;t<particles.size();t++)
         *{
         *        //printf("2\t\%lf\t\%lf\t%lf\n",
         *                                        //(*iter)->r[0],(*iter)->r[1],(*iter)->r[2]);
         *        iter++;
         *}
         */
}

void InitBarostatFull(std::list<Particle*>& particles)
{
        unsigned int i;
        unsigned int Num;
        //double sigma=1.1;
        double sigma=AmbientTemp;
        double lambda = 0;
        double vel = 0;

        Particle *px0;
        Particle *pxL;
        Particle *py0;
        Particle *pyL;
        Particle *pz0;
        Particle *pzL;

    /*
     *double vSum = 0;
     *for(int i=0;i<N;i++)
     *    for(int k=0;k<3;k++)
     *        vSum += cube[i].v[k] * cube[i].v[k];
     *double lambda = sqrt(3*(N-1)*Temp/vSum);
     *for(int i=0;i<N;i++)
     *    for(int k=0;k<3;k++)
     *        cube[i].v[k] *= lambda;
     */

        double kinE = 0.;

        //side: x0
        Num = NumberOfParticles();
        for(i=0;i<Num;i++)
        {
            
            px0 = new Particle(-L,
                            ((gsl_rng_uniform(r)-0.5)*3*L)+0.5*L,
                            ((gsl_rng_uniform(r)-0.5)*3*L)+0.5*L,
                            gsl_ran_rayleigh(r,sigma),
                            ((gsl_rng_uniform(r)-0.5)*4)*AmbientTemp,
                            ((gsl_rng_uniform(r)-0.5)*4)*AmbientTemp,"GasIn");
            /*
             *vel = px0->v[0]*px0->v[0]+px0->v[1]*px0->v[2]+px0->v[2]*px0->v[2];
             *vel = vel/3.;
             *lambda = sqrt(AmbientTemp/vel);
             *px0->v[0] = px0->v[0]*lambda;  
             *px0->v[1] = px0->v[1]*lambda;  
             *px0->v[2] = px0->v[2]*lambda;  
             */
            particles.push_back(px0);
            kinE = px0->m*(px0->v[0]*px0->v[0]+px0->v[1]*px0->v[1]+px0->v[2]*px0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side xL
        Num = NumberOfParticles();
        for(i=0;i<Num;i++)
        {
            pxL = new Particle(2.*L,
                            ((gsl_rng_uniform(r)-0.5)*3*L)+0.5*L,
                            ((gsl_rng_uniform(r)-0.5)*3*L)+0.5*L,
                            gsl_ran_rayleigh(r,sigma)*(-1.),
                            ((gsl_rng_uniform(r)-0.5)*4)*AmbientTemp,
                            ((gsl_rng_uniform(r)-0.5)*4)*AmbientTemp,"GasIn");
            /*
             *vel = pxL->v[0]*pxL->v[0]+pxL->v[1]*pxL->v[2]+pxL->v[2]*pxL->v[2];
             *vel = vel/3.;
             *lambda = sqrt(AmbientTemp/vel);
             *pxL->v[0] = pxL->v[0]*lambda;  
             *pxL->v[1] = pxL->v[1]*lambda;  
             *pxL->v[2] = pxL->v[2]*lambda;  
             */
            particles.push_back(pxL);
            kinE = pxL->m*(pxL->v[0]*pxL->v[0]+pxL->v[1]*pxL->v[1]+pxL->v[2]*pxL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }   

        //side y0
        Num = NumberOfParticles();
        for(i=0;i<Num;i++)
        {
            py0 = new Particle(((gsl_rng_uniform(r)-0.5)*3*L)+0.5*L,
                            -L,
                            ((gsl_rng_uniform(r)-0.5)*3*L)+0.5*L,
                            ((gsl_rng_uniform(r)-0.5)*4)*AmbientTemp,
                            gsl_ran_rayleigh(r,sigma),
                            ((gsl_rng_uniform(r)-0.5)*4)*AmbientTemp,"GasIn");
            /*
             *vel = py0->v[0]*py0->v[0]+py0->v[1]*py0->v[2]+py0->v[2]*py0->v[2];
             *vel = vel/3.;
             *lambda = sqrt(AmbientTemp/vel);
             *py0->v[0] = py0->v[0]*lambda;  
             *py0->v[1] = py0->v[1]*lambda;  
             *py0->v[2] = py0->v[2]*lambda;  
             */
            particles.push_back(py0);
            kinE = py0->m*(py0->v[0]*py0->v[0]+py0->v[1]*py0->v[1]+py0->v[2]*py0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }   

        //side yL
        Num = NumberOfParticles();
        for(i=0;i<Num;i++)
        {
            pyL = new Particle(((gsl_rng_uniform(r)-0.5)*3*L)+0.5*L,
                            2.*L,
                            ((gsl_rng_uniform(r)-0.5)*3*L)+0.5*L,
                            ((gsl_rng_uniform(r)-0.5)*4)*AmbientTemp,
                            gsl_ran_rayleigh(r,sigma)*(-1.),
                            ((gsl_rng_uniform(r)-0.5)*4)*AmbientTemp,"GasIn");
            /*
             *vel = pyL->v[0]*pyL->v[0]+pyL->v[1]*pyL->v[2]+pyL->v[2]*pyL->v[2];
             *vel = vel/3.;
             *lambda = sqrt(AmbientTemp/vel);
             *pyL->v[0] = pyL->v[0]*lambda;  
             *pyL->v[1] = pyL->v[1]*lambda;  
             *pyL->v[2] = pyL->v[2]*lambda;  
             */
            particles.push_back(pyL);
            kinE = pyL->m*(pyL->v[0]*pyL->v[0]+pyL->v[1]*pyL->v[1]+pyL->v[2]*pyL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }   
        
        //side z0
        Num = NumberOfParticles();
        for(i=0;i<Num;i++)
        {
            pz0 = new Particle(
                            ((gsl_rng_uniform(r)-0.5)*3*L)+0.5*L,
                            ((gsl_rng_uniform(r)-0.5)*3*L)+0.5*L,
                            -L,
                            ((gsl_rng_uniform(r)-0.5)*4)*AmbientTemp,
                            ((gsl_rng_uniform(r)-0.5)*4)*AmbientTemp,
                            gsl_ran_rayleigh(r,sigma),"GasIn");
            /*
             *vel = pz0->v[0]*pz0->v[0]+pz0->v[1]*pz0->v[2]+pz0->v[2]*pz0->v[2];
             *vel = vel/3.;
             *lambda = sqrt(AmbientTemp/vel);
             *pz0->v[0] = pz0->v[0]*lambda;  
             *pz0->v[1] = pz0->v[1]*lambda;  
             *pz0->v[2] = pz0->v[2]*lambda;  
             */
            particles.push_back(pz0);
            kinE = pz0->m*(pz0->v[0]*pz0->v[0]+pz0->v[1]*pz0->v[1]+pz0->v[2]*pz0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }   

        //side zL
        Num = NumberOfParticles();
        for(i=0;i<Num;i++)
        {
            pzL = new Particle(
                            ((gsl_rng_uniform(r)-0.5)*3*L)+0.5*L,
                            ((gsl_rng_uniform(r)-0.5)*3*L)+0.5*L,
                            2.*L,
                            ((gsl_rng_uniform(r)-0.5)*4)*AmbientTemp,
                            ((gsl_rng_uniform(r)-0.5)*4)*AmbientTemp,
                            gsl_ran_rayleigh(r,sigma)*(-1.),"GasIn");
            /*
             *vel = pzL->v[0]*pzL->v[0]+pzL->v[1]*pzL->v[2]+pzL->v[2]*pzL->v[2];
             *vel = vel/3.;
             *lambda = sqrt(AmbientTemp/vel);
             *pzL->v[0] = pzL->v[0]*lambda;  
             *pzL->v[1] = pzL->v[1]*lambda;  
             *pzL->v[2] = pzL->v[2]*lambda;  
             */
            particles.push_back(pzL);
            kinE = pzL->m*(pzL->v[0]*pzL->v[0]+pzL->v[1]*pzL->v[1]+pzL->v[2]*pzL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }   
/*
 *
 *        //side: x0
 *        N = NumberOfParticles();
 *        for(i=0;i<N;i++)
 *        {
 *            
 *            px0 = new Particle(-L/2.,
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            gsl_ran_rayleigh(r,sigma),
 *                            gsl_ran_gaussian(r,0.01)+sigma,
 *                            gsl_ran_gaussian(r,0.01)+sigma,"GasIn");
 *            particles.push_back(px0);
 *            kinE = px0->m*(px0->v[0]*px0->v[0]+px0->v[1]*px0->v[1]+px0->v[2]*px0->v[2])*0.5;
 *            gsl_histogram_increment(gas_in,kinE);
 *        }
 *
 *        //side xL
 *        N = NumberOfParticles();
 *        for(i=0;i<N;i++)
 *        {
 *            pxL = new Particle(L+L/2.,
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            gsl_ran_rayleigh(r,sigma)*(-1.),
 *                            gsl_ran_gaussian(r,0.01)+sigma,
 *                            gsl_ran_gaussian(r,0.01)+sigma,"GasIn");
 *            particles.push_back(pxL);
 *            kinE = pxL->m*(pxL->v[0]*pxL->v[0]+pxL->v[1]*pxL->v[1]+pxL->v[2]*pxL->v[2])*0.5;
 *            gsl_histogram_increment(gas_in,kinE);
 *        }   
 *
 *        //side y0
 *        N = NumberOfParticles();
 *        for(i=0;i<N;i++)
 *        {
 *            py0 = new Particle(((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            -L/2.,
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            gsl_ran_gaussian(r,0.01)+sigma,
 *                            gsl_ran_rayleigh(r,sigma),
 *                            gsl_ran_gaussian(r,0.01)+sigma,"GasIn");
 *            particles.push_back(py0);
 *            kinE = py0->m*(py0->v[0]*py0->v[0]+py0->v[1]*py0->v[1]+py0->v[2]*py0->v[2])*0.5;
 *            gsl_histogram_increment(gas_in,kinE);
 *        }   
 *
 *        //side yL
 *        N = NumberOfParticles();
 *        for(i=0;i<N;i++)
 *        {
 *            pyL = new Particle(((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            L+L/2.,
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            gsl_ran_gaussian(r,0.01)+sigma,
 *                            gsl_ran_rayleigh(r,sigma)*(-1.),
 *                            gsl_ran_gaussian(r,0.01)+sigma,"GasIn");
 *            particles.push_back(pyL);
 *            kinE = pyL->m*(pyL->v[0]*pyL->v[0]+pyL->v[1]*pyL->v[1]+pyL->v[2]*pyL->v[2])*0.5;
 *            gsl_histogram_increment(gas_in,kinE);
 *        }   
 *        
 *        //side z0
 *        N = NumberOfParticles();
 *        for(i=0;i<N;i++)
 *        {
 *            pz0 = new Particle(
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            -L/2.,
 *                            gsl_ran_gaussian(r,0.01)+sigma,
 *                            gsl_ran_gaussian(r,0.01)+sigma,
 *                            gsl_ran_rayleigh(r,sigma),"GasIn");
 *            particles.push_back(pz0);
 *            kinE = pz0->m*(pz0->v[0]*pz0->v[0]+pz0->v[1]*pz0->v[1]+pz0->v[2]*pz0->v[2])*0.5;
 *            gsl_histogram_increment(gas_in,kinE);
 *        }   
 *
 *        //side zL
 *        N = NumberOfParticles();
 *        for(i=0;i<N;i++)
 *        {
 *            pzL = new Particle(
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            L+L/2.,
 *                            gsl_ran_gaussian(r,0.01)+sigma,
 *                            gsl_ran_gaussian(r,0.01)+sigma,
 *                            gsl_ran_rayleigh(r,sigma)*(-1.),"GasIn");
 *            particles.push_back(pzL);
 *            kinE = pzL->m*(pzL->v[0]*pzL->v[0]+pzL->v[1]*pzL->v[1]+pzL->v[2]*pzL->v[2])*0.5;
 *            gsl_histogram_increment(gas_in,kinE);
 *        }   
 */
        /*
         *delete px0;
         *delete pxL;
         *delete py0;
         *delete pyL;
         *delete pz0;
         *delete pzL;
         */

        //delete px0; segfault
        //delete pxL; NO segfault
        //delete py0; NO segfault
        //delete pyL; NO segfault
        //delete pz0; segfault
        //delete pzL; NO segfault
}

void InitBarostatFullNew(std::list<Particle*>& particles)
{
    /* fix mistake that i made: forgot to add some sides
     * so for every face, there are two more "on the side"
     */

        unsigned int i;
        unsigned int N;
        double sigma=AmbientTemp;

        Particle *px0;
        Particle *pxL;
        Particle *py0;
        Particle *pyL;
        Particle *pz0;
        Particle *pzL;

        Particle *px0Left;
        Particle *px0Right;
        Particle *px0Front;
        Particle *px0Back;

        Particle *pxLLeft;
        Particle *pxLRight;
        Particle *pxLFront;
        Particle *pxLBack;

        Particle *py0Left;
        Particle *py0Right;
        Particle *py0Front;
        Particle *py0Back;

        Particle *pyLLeft;
        Particle *pyLRight;
        Particle *pyLFront;
        Particle *pyLBack;

        Particle *pz0Left;
        Particle *pz0Right;
        Particle *pz0Front;
        Particle *pz0Back;

        Particle *pzLLeft;
        Particle *pzLRight;
        Particle *pzLFront;
        Particle *pzLBack;

        double kinE = 0.;

        //side: x0
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            
            px0 = new Particle(-L*0.5,
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(px0);
            //kinE = px0->m*(px0->v[0]*px0->v[0]+px0->v[1]*px0->v[1]+px0->v[2]*px0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side xL
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pxL = new Particle(L+0.5*L,
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma)*(-1.),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pxL);
            //kinE = pxL->m*(pxL->v[0]*pxL->v[0]+pxL->v[1]*pxL->v[1]+pxL->v[2]*pxL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }   

        //side y0
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            py0 = new Particle(gsl_rng_uniform(r)*L,
                            -L*0.5,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(py0);
            //kinE = py0->m*(py0->v[0]*py0->v[0]+py0->v[1]*py0->v[1]+py0->v[2]*py0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }   

        //side yL
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pyL = new Particle(gsl_rng_uniform(r)*L,
                            L+0.5*L,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma)*(-1.),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pyL);
            //kinE = pyL->m*(pyL->v[0]*pyL->v[0]+pyL->v[1]*pyL->v[1]+pyL->v[2]*pyL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }   
        
        //side z0
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pz0 = new Particle(
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            -L*0.5,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma),"GasIn");
            particles.push_back(pz0);
            //kinE = pz0->m*(pz0->v[0]*pz0->v[0]+pz0->v[1]*pz0->v[1]+pz0->v[2]*pz0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }   

        //side zL
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pzL = new Particle(
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            L+0.5*L,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma)*(-1.),"GasIn");
            particles.push_back(pzL);
            //kinE = pzL->m*(pzL->v[0]*pzL->v[0]+pzL->v[1]*pzL->v[1]+pzL->v[2]*pzL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }   

        //side: x0Left
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            px0Left = new Particle((gsl_rng_uniform(r)*(-1.)*0.5*L),
                            L,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(px0Left);
            //kinE = px0Left->m*(px0->v[0]*px0->v[0]+px0->v[1]*px0->v[1]+px0->v[2]*px0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: x0Right
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            px0Right = new Particle((gsl_rng_uniform(r)*(-1.)*0.5*L),
                            0.,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(px0Right);
            //kinE = px0Right->m*(px0->v[0]*px0->v[0]+px0->v[1]*px0->v[1]+px0->v[2]*px0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: x0Back
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            px0Back = new Particle((gsl_rng_uniform(r)*(-1.)*0.5*L),
                            gsl_rng_uniform(r)*L,
                            0.,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(px0Back);
            //kinE = px0Back->m*(px0->v[0]*px0->v[0]+px0->v[1]*px0->v[1]+px0->v[2]*px0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: x0Front
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            px0Front = new Particle((gsl_rng_uniform(r)*(-1.)*0.5*L),
                            gsl_rng_uniform(r)*L,
                            L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(px0Front);
            //kinE = px0Front->m*(px0->v[0]*px0->v[0]+px0->v[1]*px0->v[1]+px0->v[2]*px0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: xLLeft
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pxLLeft = new Particle((gsl_rng_uniform(r)*0.5*L)+L,
                            L,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pxLLeft);
            //kinE = pxLLeft->m*(pxL->v[0]*pxL->v[0]+pxL->v[1]*pxL->v[1]+pxL->v[2]*pxL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: xLRight
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pxLRight = new Particle((gsl_rng_uniform(r)*0.5*L)+L,
                            0.,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pxLRight);
            //kinE = pxLRight->m*(pxL->v[0]*pxL->v[0]+pxL->v[1]*pxL->v[1]+pxL->v[2]*pxL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: xLBack
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pxLBack = new Particle((gsl_rng_uniform(r)*0.5*L)+L,
                            gsl_rng_uniform(r)*L,
                            0.,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pxLBack);
            //kinE = pxLBack->m*(pxL->v[0]*pxL->v[0]+pxL->v[1]*pxL->v[1]+pxL->v[2]*pxL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: xLFront
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pxLFront = new Particle((gsl_rng_uniform(r)*0.5*L)+L,
                            gsl_rng_uniform(r)*L,
                            L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pxLFront);
            //kinE = pxLFront->m*(pxL->v[0]*pxL->v[0]+pxL->v[1]*pxL->v[1]+pxL->v[2]*pxL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: y0Left
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            py0Left = new Particle(0.,
                            (gsl_rng_uniform(r)*(-1.)*0.5*L),
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(py0Left);
            //kinE = py0Left->m*(py0->v[0]*py0->v[0]+py0->v[1]*py0->v[1]+py0->v[2]*py0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: y0Right
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            py0Right = new Particle(L,
                            (gsl_rng_uniform(r)*(-1.)*0.5*L),
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(py0Right);
            //kinE = py0Right->m*(py0->v[0]*py0->v[0]+py0->v[1]*py0->v[1]+py0->v[2]*py0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }
        
        //side: y0Back
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            py0Back = new Particle(gsl_rng_uniform(r)*L,
                            (gsl_rng_uniform(r)*(-1.)*0.5*L),
                            0.,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(py0Back);
            //kinE = py0Back->m*(py0->v[0]*py0->v[0]+py0->v[1]*py0->v[1]+py0->v[2]*py0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: y0Front
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            py0Front = new Particle(gsl_rng_uniform(r)*L,
                            (gsl_rng_uniform(r)*(-1.)*0.5*L),
                            L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(py0Front);
            //kinE = py0Front->m*(py0->v[0]*py0->v[0]+py0->v[1]*py0->v[1]+py0->v[2]*py0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }
        
        //side: yLLeft
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pyLLeft = new Particle(0.,
                            (gsl_rng_uniform(r)*0.5*L)+L,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pyLLeft);
            //kinE = pyLLeft->m*(pyL->v[0]*pyL->v[0]+pyL->v[1]*pyL->v[1]+pyL->v[2]*pyL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: yLRight
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pyLRight = new Particle(L,
                            (gsl_rng_uniform(r)*0.5*L)+L,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pyLRight);
            //kinE = pyLRight->m*(pyL->v[0]*pyL->v[0]+pyL->v[1]*pyL->v[1]+pyL->v[2]*pyL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }
        
        //side: yLBack
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pyLBack = new Particle(gsl_rng_uniform(r)*L,
                            (gsl_rng_uniform(r)*0.5*L)+L,
                            0.,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pyLBack);
            //kinE = pyLBack->m*(pyL->v[0]*pyL->v[0]+pyL->v[1]*pyL->v[1]+pyL->v[2]*pyL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: yLFront
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pyLFront = new Particle(gsl_rng_uniform(r)*L,
                            (gsl_rng_uniform(r)*0.5*L)+L,
                            L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pyLFront);
            //kinE = pyLFront->m*(pyL->v[0]*pyL->v[0]+pyL->v[1]*pyL->v[1]+pyL->v[2]*pyL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: z0Left
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pz0Left = new Particle(0.,
                            gsl_rng_uniform(r)*L,
                            (gsl_rng_uniform(r)*(-1.)*0.5*L),
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pz0Left);
            //kinE = pz0Left->m*(pz0->v[0]*pz0->v[0]+pz0->v[1]*pz0->v[1]+pz0->v[2]*pz0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: z0Right
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pz0Right = new Particle(L,
                            gsl_rng_uniform(r)*L,
                            (gsl_rng_uniform(r)*(-1.)*0.5*L),
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pz0Right);
            //kinE = pz0Right->m*(pz0->v[0]*pz0->v[0]+pz0->v[1]*pz0->v[1]+pz0->v[2]*pz0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }
        
        //side: z0Back
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pz0Back = new Particle(gsl_rng_uniform(r)*L,
                            0.,
                            (gsl_rng_uniform(r)*(-1.)*0.5*L),
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pz0Back);
            //kinE = pz0Back->m*(pz0->v[0]*pz0->v[0]+pz0->v[1]*pz0->v[1]+pz0->v[2]*pz0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: z0Front
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pz0Front = new Particle(gsl_rng_uniform(r)*L,
                            L,
                            (gsl_rng_uniform(r)*(-1.)*0.5*L),
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pz0Front);
            //kinE = pz0Front->m*(pz0->v[0]*pz0->v[0]+pz0->v[1]*pz0->v[1]+pz0->v[2]*pz0->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }
        
        //side: zLLeft
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pzLLeft = new Particle(0.,
                            gsl_rng_uniform(r)*L,
                            (gsl_rng_uniform(r)*0.5*L)+L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pzLLeft);
            //kinE = pzLLeft->m*(pzL->v[0]*pzL->v[0]+pzL->v[1]*pzL->v[1]+pzL->v[2]*pzL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: zLRight
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pzLRight = new Particle(L,
                            gsl_rng_uniform(r)*L,
                            (gsl_rng_uniform(r)*0.5*L)+L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pzLRight);
            //kinE = pzLRight->m*(pzL->v[0]*pzL->v[0]+pzL->v[1]*pzL->v[1]+pzL->v[2]*pzL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }
        
        //side: zLBack
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pzLBack = new Particle(gsl_rng_uniform(r)*L,
                            0.,
                            (gsl_rng_uniform(r)*0.5*L)+L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pzLBack);
            //kinE = pzLBack->m*(pzL->v[0]*pzL->v[0]+pzL->v[1]*pzL->v[1]+pzL->v[2]*pzL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }

        //side: zLFront
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pzLFront = new Particle(gsl_rng_uniform(r)*L,
                            L,
                            (gsl_rng_uniform(r)*0.5*L)+L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pzLFront);
            //kinE = pzLFront->m*(pzL->v[0]*pzL->v[0]+pzL->v[1]*pzL->v[1]+pzL->v[2]*pzL->v[2])*0.5;
            //gsl_histogram_increment(gas_in,kinE);
        }
/*
 *
 *        //side: x0
 *        N = NumberOfParticles();
 *        for(i=0;i<N;i++)
 *        {
 *            
 *            px0 = new Particle(-L/2.,
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            gsl_ran_rayleigh(r,sigma),
 *                            gsl_ran_gaussian(r,sigma),
 *                            gsl_ran_gaussian(r,sigma),"GasIn");
 *            particles.push_back(px0);
 *            kinE = px0->m*(px0->v[0]*px0->v[0]+px0->v[1]*px0->v[1]+px0->v[2]*px0->v[2])*0.5;
 *            gsl_histogram_increment(gas_in,kinE);
 *        }
 *
 *        //side xL
 *        N = NumberOfParticles();
 *        for(i=0;i<N;i++)
 *        {
 *            pxL = new Particle(L+L/2.,
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            gsl_ran_rayleigh(r,sigma)*(-1.),
 *                            gsl_ran_gaussian(r,sigma),
 *                            gsl_ran_gaussian(r,sigma),"GasIn");
 *            particles.push_back(pxL);
 *            kinE = pxL->m*(pxL->v[0]*pxL->v[0]+pxL->v[1]*pxL->v[1]+pxL->v[2]*pxL->v[2])*0.5;
 *            gsl_histogram_increment(gas_in,kinE);
 *        }   
 *
 *        //side y0
 *        N = NumberOfParticles();
 *        for(i=0;i<N;i++)
 *        {
 *            py0 = new Particle(((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            -L/2.,
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            gsl_ran_gaussian(r,sigma),
 *                            gsl_ran_rayleigh(r,sigma),
 *                            gsl_ran_gaussian(r,sigma),"GasIn");
 *            particles.push_back(py0);
 *            kinE = py0->m*(py0->v[0]*py0->v[0]+py0->v[1]*py0->v[1]+py0->v[2]*py0->v[2])*0.5;
 *            gsl_histogram_increment(gas_in,kinE);
 *        }   
 *
 *        //side yL
 *        N = NumberOfParticles();
 *        for(i=0;i<N;i++)
 *        {
 *            pyL = new Particle(((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            L+L/2.,
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            gsl_ran_gaussian(r,sigma),
 *                            gsl_ran_rayleigh(r,sigma)*(-1.),
 *                            gsl_ran_gaussian(r,sigma),"GasIn");
 *            particles.push_back(pyL);
 *            kinE = pyL->m*(pyL->v[0]*pyL->v[0]+pyL->v[1]*pyL->v[1]+pyL->v[2]*pyL->v[2])*0.5;
 *            gsl_histogram_increment(gas_in,kinE);
 *        }   
 *        
 *        //side z0
 *        N = NumberOfParticles();
 *        for(i=0;i<N;i++)
 *        {
 *            pz0 = new Particle(
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            -L/2.,
 *                            gsl_ran_gaussian(r,sigma),
 *                            gsl_ran_gaussian(r,sigma),
 *                            gsl_ran_rayleigh(r,sigma),"GasIn");
 *            particles.push_back(pz0);
 *            kinE = pz0->m*(pz0->v[0]*pz0->v[0]+pz0->v[1]*pz0->v[1]+pz0->v[2]*pz0->v[2])*0.5;
 *            gsl_histogram_increment(gas_in,kinE);
 *        }   
 *
 *        //side zL
 *        N = NumberOfParticles();
 *        for(i=0;i<N;i++)
 *        {
 *            pzL = new Particle(
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            ((gsl_rng_uniform(r)-0.5)*2*L)+0.5*L,
 *                            L+L/2.,
 *                            gsl_ran_gaussian(r,sigma),
 *                            gsl_ran_gaussian(r,sigma),
 *                            gsl_ran_rayleigh(r,sigma)*(-1.),"GasIn");
 *            particles.push_back(pzL);
 *            kinE = pzL->m*(pzL->v[0]*pzL->v[0]+pzL->v[1]*pzL->v[1]+pzL->v[2]*pzL->v[2])*0.5;
 *            gsl_histogram_increment(gas_in,kinE);
 *        }   
 */
        /*
         *delete px0;
         *delete pxL;
         *delete py0;
         *delete pyL;
         *delete pz0;
         *delete pzL;
         */

        //delete px0; segfault
        //delete pxL; NO segfault
        //delete py0; NO segfault
        //delete pyL; NO segfault
        //delete pz0; segfault
        //delete pzL; NO segfault
}

void InitBarostatNew(std::list<Particle*>& particles)
{
        //FILE *baroOut = fopen("output/barostat.dat","a");
        //int i,j,k;
        unsigned int i;
        unsigned int N;
        double sigma=1.1;

        Particle *px0;
        Particle *pxL;
        Particle *py0;
        Particle *pyL;
        Particle *pz0;
        Particle *pzL;

        double kinE = 0.;

        //side: x0
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            
            px0 = new Particle(-L/2.,
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(px0);
            kinE = px0->m*(px0->v[0]*px0->v[0]+px0->v[1]*px0->v[1]+px0->v[2]*px0->v[2])*0.5;
            gsl_histogram_increment(gas_in,kinE);
        }

        //side xL
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pxL = new Particle(L+L/2.,
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma)*(-1.),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pxL);
            kinE = pxL->m*(pxL->v[0]*pxL->v[0]+pxL->v[1]*pxL->v[1]+pxL->v[2]*pxL->v[2])*0.5;
            gsl_histogram_increment(gas_in,kinE);
        }   

        //side y0
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            py0 = new Particle(gsl_rng_uniform(r)*L,
                            -L/2.,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(py0);
            kinE = py0->m*(py0->v[0]*py0->v[0]+py0->v[1]*py0->v[1]+py0->v[2]*py0->v[2])*0.5;
            gsl_histogram_increment(gas_in,kinE);
        }   

        //side yL
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pyL = new Particle(gsl_rng_uniform(r)*L,
                            L+L/2.,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma)*(-1.),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pyL);
            kinE = pyL->m*(pyL->v[0]*pyL->v[0]+pyL->v[1]*pyL->v[1]+pyL->v[2]*pyL->v[2])*0.5;
            gsl_histogram_increment(gas_in,kinE);
        }   
        
        //side z0
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pz0 = new Particle(
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            -L/2.,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma),"GasIn");
            particles.push_back(pz0);
            kinE = pz0->m*(pz0->v[0]*pz0->v[0]+pz0->v[1]*pz0->v[1]+pz0->v[2]*pz0->v[2])*0.5;
            gsl_histogram_increment(gas_in,kinE);
        }   

        //side zL
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pzL = new Particle(
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            L+L/2.,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma)*(-1.),"GasIn");
            particles.push_back(pzL);
            kinE = pzL->m*(pzL->v[0]*pzL->v[0]+pzL->v[1]*pzL->v[1]+pzL->v[2]*pzL->v[2])*0.5;
            gsl_histogram_increment(gas_in,kinE);
        }   
        /*
         *delete px0;
         *delete pxL;
         *delete py0;
         *delete pyL;
         *delete pz0;
         *delete pzL;
         */

        //delete px0; segfault
        //delete pxL; NO segfault
        //delete py0; NO segfault
        //delete pyL; NO segfault
        //delete pz0; segfault
        //delete pzL; NO segfault
        
}

void InitBarostat(std::list<Particle*>& particles,std::list<Particle*>& particlesHistory)
{
        //FILE *baroOut = fopen("output/barostat.dat","a");
        //int i,j,k;
        unsigned int i;
        unsigned int N;
        double sigma=1.1;

        Particle *px0;
        Particle *pxL;
        Particle *py0;
        Particle *pyL;
        Particle *pz0;
        Particle *pzL;

        //side: x0
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            
            px0 = new Particle(-L/2.,
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(px0);
            particlesHistory.push_back(px0);
        }

        //side xL
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pxL = new Particle(L+L/2.,
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_rayleigh(r,sigma)*(-1.),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pxL);
            particlesHistory.push_back(pxL);
        }   

        //side y0
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            py0 = new Particle(gsl_rng_uniform(r)*L,
                            -L/2.,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(py0);
            particlesHistory.push_back(py0);
        }   

        //side yL
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pyL = new Particle(gsl_rng_uniform(r)*L,
                            L+L/2.,
                            gsl_rng_uniform(r)*L,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma)*(-1.),
                            gsl_ran_gaussian(r,sigma),"GasIn");
            particles.push_back(pyL);
            particlesHistory.push_back(pyL);
        }   
        
        //side z0
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pz0 = new Particle(
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            -L/2.,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma),"GasIn");
            particles.push_back(pz0);
            particlesHistory.push_back(pz0);
        }   

        //side zL
        N = NumberOfParticles();
        for(i=0;i<N;i++)
        {
            pzL = new Particle(
                            gsl_rng_uniform(r)*L,
                            gsl_rng_uniform(r)*L,
                            L+L/2.,
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_gaussian(r,sigma),
                            gsl_ran_rayleigh(r,sigma)*(-1.),"GasIn");
            particles.push_back(pzL);
            particlesHistory.push_back(pzL);
        }   
        /*
         *delete px0;
         *delete pxL;
         *delete py0;
         *delete pyL;
         *delete pz0;
         *delete pzL;
         */

        //delete px0; segfault
        //delete pxL; NO segfault
        //delete py0; NO segfault
        //delete pyL; NO segfault
        //delete pz0; segfault
        //delete pzL; NO segfault
        
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

void InitBarostat(std::vector<Particle*>& particles,std::vector<Particle*>& particleHistory)
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

        //std::string filename;
        //std::string filename = "output/combined/combined_";
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
        fprintf(combinedOut,"%ld\nFrame 0\n",N+gas.size());
        std::list<Particle*>::iterator gasIter = gas.begin();
        unsigned int t=0;
        for(unsigned int j=0;j<N;j++)
        {
                fprintf(combinedOut,"%s\t\%lf\t\%lf\t%lf\t%d\n",(cube[j].name).c_str(),cube[j].r[0],cube[j].r[1],cube[j].r[2],cube[j].ID);
                //printed++;
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
        //printDistances(cube,gas,Run);
        
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
    //double lowerBound = -(L + eps);
    //double upperBound = 2.*L + eps;

    for(temp=particles.begin();temp!=particles.end();)
    {
        if((*temp)->r[0] < lowerBound || (*temp)->r[0] > upperBound || 
                (*temp)->r[1] < lowerBound || (*temp)->r[1] > upperBound 
                ||(*temp)->r[2] < lowerBound || (*temp)->r[2] > upperBound)
        {
            //std::cout << "Deleted [ID]: " << (*temp)->ID << std::endl;
            temp=particles.erase(temp);
            //(*temp)->name = "GasGone";
        }
        else
            ++temp;
    }

}

void CheckBoundariesNoCube(std::list<Particle*> &particles)
{
    std::list<Particle*>::iterator temp;
    double lowerBound = -(L/2 + eps);
    double upperBound = L + L/2 + eps;
    //double lowerBound = -(L + eps);
    //double upperBound = 2.*L + eps;

    for(temp=particles.begin();temp!=particles.end();)
    {
        if((*temp)->r[0] < lowerBound || (*temp)->r[0] > upperBound || 
                (*temp)->r[1] < lowerBound || (*temp)->r[1] > upperBound 
                ||(*temp)->r[2] < lowerBound || (*temp)->r[2] > upperBound)
        {
            //std::cout << "Deleted [ID]: " << (*temp)->ID << std::endl;
            temp=particles.erase(temp);
            //(*temp)->name = "GasGone";
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
    double aTotOld = 0;
    double aTotNew = 0;
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
        aTotOld = sqrt((*temp)->aOld[0]*(*temp)->aOld[0]+(*temp)->aOld[1]*(*temp)->aOld[1]+(*temp)->aOld[2]*(*temp)->aOld[2]);
        aTotNew = sqrt((*temp)->a[0]*(*temp)->a[0]+(*temp)->a[1]*(*temp)->a[1]+(*temp)->a[2]*(*temp)->a[2]);
        if(aTotOld > aTotNew)
            (*temp)->name = "GasOut";
        if((*temp)->r[0] < lowerBound || (*temp)->r[0] > upperBound || 
                (*temp)->r[1] < lowerBound || (*temp)->r[1] > upperBound 
                ||(*temp)->r[2] < lowerBound || (*temp)->r[2] > upperBound)
        {
            //std::cout << "Deleted [ID]: " << (*temp)->ID << std::endl;
            temp=particles.erase(temp);
            //(*temp)->name = "GasGone";
        }
        else
            ++temp;
    }

}

void CheckBoundariesNew(Particle* cube,std::list<Particle*> &particles)
{
    std::list<Particle*>::iterator temp;
    //double lowerBound = -(L/2 + eps);
    //double upperBound = L + L/2 + eps;
    double lowerBound = -(L + eps);
    double upperBound = 2.*L + eps;
    double aTotOld = 0;
    double aTotNew = 0;
    //double dist[3] = {0,0,0};
    //double distSqd = 0.;
    //double vSqd = 0.;

    for(temp=particles.begin();temp!=particles.end();)
    {
        double kinE = 0.;
        double kinEOld = 0;
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
        aTotOld = sqrt((*temp)->aOld[0]*(*temp)->aOld[0]+(*temp)->aOld[1]*(*temp)->aOld[1]+(*temp)->aOld[2]*(*temp)->aOld[2]);
        aTotNew = sqrt((*temp)->a[0]*(*temp)->a[0]+(*temp)->a[1]*(*temp)->a[1]+(*temp)->a[2]*(*temp)->a[2]);

        /*
         *if(aTotOld > aTotNew)
         *    (*temp)->name = "GasOut";
         */
        if(aTotNew > aTotOld && strcmp((*temp)->name.c_str(),"GasIn") == 0)
        {
            (*temp)->name = "GasChange";
            kinEOld = (*temp)->m*((*temp)->vOld[0]*(*temp)->vOld[0]+(*temp)->vOld[1]*(*temp)->vOld[1]+(*temp)->vOld[2]*(*temp)->vOld[2])*0.5;
            gsl_histogram_increment(gas_real_in,kinEOld);
        }
        if(std::abs(aTotNew-aTotOld) < 10e-9 && strcmp((*temp)->name.c_str(),"GasChange") == 0)
        {
            (*temp)->name = "GasOut";
            kinE = (*temp)->m*((*temp)->v[0]*(*temp)->v[0]+(*temp)->v[1]*(*temp)->v[1]+(*temp)->v[2]*(*temp)->v[2])*0.5;
            gsl_histogram_increment(gas_out,kinE);
        }
        if((*temp)->r[0] < lowerBound || (*temp)->r[0] > upperBound || 
                (*temp)->r[1] < lowerBound || (*temp)->r[1] > upperBound 
                ||(*temp)->r[2] < lowerBound || (*temp)->r[2] > upperBound)
        {
            //std::cout << "Deleted [ID]: " << (*temp)->ID << std::endl;
            temp=particles.erase(temp);
            //(*temp)->name = "GasGone";
            //++temp;
        }
        else
            ++temp;
    }

}

void CheckBoundariesNoRemove(Particle* cube,std::list<Particle*> &particles)
{
    std::list<Particle*>::iterator temp;
    //double lowerBound = -(L/2 + eps);
    //double upperBound = L + L/2 + eps;
    double lowerBound = -(L + eps);
    double upperBound = 2.*L + eps;
    double aTotOld = 0;
    double aTotNew = 0;
    //double dist[3] = {0,0,0};
    //double distSqd = 0.;
    //double vSqd = 0.;

    for(temp=particles.begin();temp!=particles.end();)
    {
        double kinE = 0.;
        double kinEOld = 0;
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
        aTotOld = sqrt((*temp)->aOld[0]*(*temp)->aOld[0]+(*temp)->aOld[1]*(*temp)->aOld[1]+(*temp)->aOld[2]*(*temp)->aOld[2]);
        aTotNew = sqrt((*temp)->a[0]*(*temp)->a[0]+(*temp)->a[1]*(*temp)->a[1]+(*temp)->a[2]*(*temp)->a[2]);

        /*
         *if(aTotOld > aTotNew)
         *    (*temp)->name = "GasOut";
         */
        if(aTotNew > aTotOld && strcmp((*temp)->name.c_str(),"GasIn") == 0)
        {
            (*temp)->name = "GasChange";
            kinEOld = (*temp)->m*((*temp)->vOld[0]*(*temp)->vOld[0]+(*temp)->vOld[1]*(*temp)->vOld[1]+(*temp)->vOld[2]*(*temp)->vOld[2])*0.5;
            gsl_histogram_increment(gas_real_in,kinEOld);
        }
        if(std::abs(aTotNew-aTotOld) < 10e-9 && strcmp((*temp)->name.c_str(),"GasChange") == 0)
        {
            (*temp)->name = "GasOut";
            kinE = (*temp)->m*((*temp)->v[0]*(*temp)->v[0]+(*temp)->v[1]*(*temp)->v[1]+(*temp)->v[2]*(*temp)->v[2])*0.5;
            gsl_histogram_increment(gas_out,kinE);
        }
        temp++;
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

void Barostat(Particle* cube, std::list<Particle*>& gas,std::list<Particle*>& gasHistory)
{
    std::list<Particle*>::iterator gasIter;
    InitBarostat(gas,gasHistory);
    
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

void BarostatNew(Particle* cube, std::list<Particle*>& gas)
{
    std::list<Particle*>::iterator gasIter;
    //InitBarostatNew(gas);
    InitBarostatFull(gas);
    
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
    CheckBoundariesNew(cube,gas);
}

void BarostatNoCube(std::list<Particle*>& gas)
{
    std::list<Particle*>::iterator gasIter;
    //InitBarostatNew(gas);
    InitBarostatFull(gas);
    
    for(gasIter = gas.begin();gasIter != gas.end(); gasIter++)
    {
        for(unsigned int m=0;m<3;m++)
        {
            (*gasIter)->r[m] += (*gasIter)->v[m]*dt + 0.5 * (*gasIter)->a[m] * dt * dt;
            (*gasIter)->v[m] += 0.5 * (*gasIter)->a[m] * dt;
        }
    }

    ComputeSoftSphereNoCube(gas);
    for(gasIter = gas.begin();gasIter != gas.end(); gasIter++)
        for(unsigned int m=0;m<3;m++)
            (*gasIter)->v[m] += 0.5 * (*gasIter)->a[m] * dt;
    CheckBoundariesNoCube(gas);
}
void BarostatNoBoundaries(Particle* cube, std::list<Particle*>& gas)
{
    std::list<Particle*>::iterator gasIter;
    //InitBarostatNew(gas);
    
    for(gasIter = gas.begin();gasIter != gas.end(); gasIter++)
    {
        for(unsigned int m=0;m<3;m++)
        {
            (*gasIter)->r[m] += (*gasIter)->v[m]*dt + 0.5 * (*gasIter)->a[m] * dt * dt/(*gasIter)->m;
            (*gasIter)->v[m] += 0.5 * (*gasIter)->a[m] * dt/(*gasIter)->m;
        }
    }

    ComputeSoftSphereTest(gas,cube);
    for(gasIter = gas.begin();gasIter != gas.end(); gasIter++)
        for(unsigned int m=0;m<3;m++)
            (*gasIter)->v[m] += 0.5 * (*gasIter)->a[m] * dt/(*gasIter)->m;
    CheckBoundariesNoRemove(cube,gas);
}

void BarostatTest(Particle* cube, std::list<Particle*>& gas)
{
    std::list<Particle*>::iterator gasIter;
    //InitBarostatFullNew(gas);
    
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
    CheckBoundariesNew(cube,gas);
}

void BarostatNew(Particle* cube, std::list<Particle*>& gas,std::list<Particle*>& gasHistory)
{
    std::list<Particle*>::iterator gasIter;
    InitBarostat(gas,gasHistory);
    
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
    CheckBoundariesNew(cube,gas);
}

void ComputeSoftSphere(std::list<Particle*>& gas, Particle* cube)
{
    std::list<Particle*>::iterator iter;
    double rij[3];
    double rSqd=0;
    double f;
    /*
     *for(unsigned int i=0;i<N;i++)
     *    for(unsigned int m=0;m<3;m++)
     *        cube[i].a[m] = 0;
     */
    for (iter=gas.begin();iter!=gas.end();iter++)         // set all accelerations to zero
    {
        for(unsigned int i=0;i<3;i++)
        {
            (*iter)->aOld[i] = (*iter)->a[i];
            (*iter)->vOld[i] = (*iter)->v[i];
        }
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
                    (*iter) -> a[m] += rij[m]*f/(*iter)->m;
                    cube[n].a[m] -= rij[m] *f;
                }
            }
        }
    }
}

void ComputeSoftSphereNoCube(std::list<Particle*>& gas)
{
    std::list<Particle*>::iterator iter;
    double rij[3];
    double rSqd=0;
    double f;
    /*
     *for(unsigned int i=0;i<N;i++)
     *    for(unsigned int m=0;m<3;m++)
     *        cube[i].a[m] = 0;
     */
    for (iter=gas.begin();iter!=gas.end();iter++)         // set all accelerations to zero
    {
        for(unsigned int i=0;i<3;i++)
        {
            (*iter)->aOld[i] = (*iter)->a[i];
            (*iter)->vOld[i] = (*iter)->v[i];
        }
        for(unsigned int i=0;i<3;i++)
        {
            (*iter)->a[i]=0;
            (*iter)->type = 3;
        }
    }
}

void ComputeSoftSphereTest(std::list<Particle*>& gas, Particle* cube)
{
    std::list<Particle*>::iterator iter;
    double rij[3];
    double rSqd=0;
    double f;
    /*
     *for(unsigned int i=0;i<N;i++)
     *    for(unsigned int m=0;m<3;m++)
     *        cube[i].a[m] = 0;
     */
    for (iter=gas.begin();iter!=gas.end();iter++)         // set all accelerations to zero
    {
        for(unsigned int i=0;i<3;i++)
        {
            (*iter)->aOld[i] = (*iter)->a[i];
            (*iter)->vOld[i] = (*iter)->v[i];
        }
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

//void eHEX(Particle* cube,FILE* output)
void eHEX(Particle* cube)
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
	K=K/2.;
	if(K==0)
		printf("Error, K must not be 0!");
	xi = sqrt(1+dQ/K);

	for(n=0;n<N;n++)
		for(i=0;i<3;i++)
			eta[n][i]=(FG/(2.*K))*cube[n].v[i];

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

void eHEXNew(Particle* cube)
{
	double eps[N][3];
	double K=0;
	double eta[N][3];
	double xi=0;
	double sumforces[3]={0,0,0};
	double sumvel=0; 
	unsigned int n,i; 
    //double v_up,vhalf_up,vnew_up;
    //double rnew_up;

    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            K += (cube[n].v[i]*cube[n].v[i]);
    K = K/2.;
    xi = sqrt(1+dQ/K);




    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            cube[n].vnew[i] = xi*cube[n].v[i]; //(19a)
    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            cube[n].vhalf[i] = cube[n].vnew[i]+dt*0.5*cube[n].a[i]; //(19b)
    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            cube[n].rnew[i] = cube[n].r[i]+dt*cube[n].vhalf[i]; //(19c)
    for(n=0;n<N;n++)
    {
        for(i=0;i<3;i++)
        {
            cube[n].v[i] = cube[n].vhalf[i];
            cube[n].r[i] = cube[n].rnew[i];
        }
    }
    ComputeAccelerations(cube); //(19d)
    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            cube[n].vnew[i] = cube[n].vhalf[i] + dt*0.5*cube[n].a[i]; //(19e)
    K = 0; 
    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            K += (cube[n].vnew[i]*cube[n].vnew[i]);
    K = K/2.;
    xi = sqrt(1+dQ/K); //\bar{xi} from (19f)
    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            cube[n].v[i] = xi*cube[n].vnew[i]; //(19f)

	for(n=0;n<N;n++)
		for(i=0;i<3;i++)
			eta[n][i]=(FG/(2.*K))*cube[n].vnew[i];
    sumvel=0;
    for(i=0;i<3;i++)
        sumforces[i]=0;
	for(n=0;n<N;n++)
		for(i=0;i<3;i++)
		{
			sumforces[i]+=cube[n].a[i];
			sumvel+=cube[n].vnew[i]*cube[n].a[i];
		}
	for(n=0;n<N;n++)
		for(i=0;i<3;i++)
            eps[n][i]=(eta[n][i]/K)*((FG/48.)+(1./6.)*sumvel)-(FG/(12.*K))*
            (cube[n].a[i]-(sumforces[i])); //(20)
    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            cube[n].r[i] = cube[n].rnew[i] - dt*dt*dt*eps[n][i]; //(19g)
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

void calcTemp(Particle* cube,FILE* output,int run)
{
    double T=0;
    for(unsigned int i=0;i<N;i++)
        for(unsigned int k=0;k<3;k++)
            T+=cube[i].v[k]*cube[i].v[k];
    T = T/(3*N);
    fprintf(output,"%d\t%lf\n",run,T);
}

void calcCOMTemp(double* vCOM,FILE* output)
{
    double T=0;
    for(unsigned int k=0;k<3;k++)
        T+=vCOM[k]*vCOM[k];
    T = T/(3);
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

void calcCM(Particle* particles,double *rCM, FILE* output)
{
    for(unsigned int i=0;i<3;i++)
        rCM[i] = 0;

    for(unsigned int i=0;i<N;i++)
        for(unsigned int k=0;k<3;k++)
            rCM[k]+=particles[i].r[k];

    for(unsigned int i=0;i<3;i++)
        rCM[i]/=N;
    fprintf(output,"%lf\t%lf\t%lf\n",rCM[0],rCM[1],rCM[2]);
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

void writePositions(Particle* particles, std::list<Particle*>gas, std::string filename)
{
    FILE* writeFile = fopen(filename.c_str(),"wb");
    std::list<Particle*>::iterator gasIter;
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
    for(gasIter=gas.begin();gasIter!=gas.end();gasIter++)
    {
        for(unsigned int j=0;j<3;j++)
            fwrite(&((*gasIter)->r[j]),sizeof(double),1,writeFile);
        for(unsigned int j=0;j<3;j++)
            fwrite(&((*gasIter)->v[j]),sizeof(double),1,writeFile);
        for(unsigned int j=0;j<3;j++)
            fwrite(&((*gasIter)->a[j]),sizeof(double),1,writeFile);
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
    const double k = 150.;

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
             fprintf(output,"%d\ntest\n",N+1);
             /*
              *fprintf(output,"%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
              *        ((*gasIter)->name).c_str(),(*gasIter)->r[0],(*gasIter)->r[1],(*gasIter)->r[2],\
              *        (*gasIter)->v[0],(*gasIter)->v[1],(*gasIter)->v[2],sqrt((*gasIter)->v[0]*(*gasIter)->v[0]+(*gasIter)->v[1]*(*gasIter)->v[1]+(*gasIter)->v[2]*(*gasIter)->v[2]),\
              *        (*gasIter)->a[0],(*gasIter)->a[1],(*gasIter)->a[2],sqrt((*gasIter)->a[0]*(*gasIter)->a[0]+(*gasIter)->a[1]*(*gasIter)->a[1]+(*gasIter)->a[2]*(*gasIter)->a[2]),\
              *        (*gasIter)->potE),(*gasIter)->m*((*gasIter)->v[0]*(*gasIter)->v[0]+(*gasIter)->v[1]*(*gasIter)->v[1]+(*gasIter)->v[2]*(*gasIter)->v[2])/2.;
              */

             /*
              * Row     variable name
              * ---------------------------
              *  1      Name
              *  2      x
              *  3      y
              *  4      z
              *  5      v_x
              *  6      v_y
              *  7      v_z
              *  8      |v|
              *  9      a_x
              *  10     a_y
              *  11     a_z      
              *  12     |a|
              *  13     E_pot
              *  14     E_tot
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
             for(unsigned int j=0;j<N;j++)
             {
                     fprintf(output,"%s\t\%lf\t\%lf\t%lf\n",(cube[j].name).c_str(),cube[j].r[0],cube[j].r[1],cube[j].r[2]);
             }


         }
     }
}

std::string DateToString()
{
    time_t now = time(0);
    tm* ltm = localtime(&now);
    std::string datestring = "";
    datestring += numberToString((1900+ltm->tm_year)-2000);
    if(1+ltm->tm_mon < 10 || 1+ltm->tm_mon == 9)
        datestring += "0";
    datestring += numberToString(1+ltm->tm_mon);
    if(ltm->tm_mday < 10)
        datestring += "0";
    datestring += numberToString(ltm->tm_mday);
    datestring += "_";
    if(ltm->tm_hour < 10)
        datestring += "0";
    datestring += numberToString(ltm->tm_hour);
    if(ltm->tm_min < 10)
        datestring += "0";
    datestring += numberToString(ltm->tm_min);
    return datestring;
     //FILE* runData = fopen(filename.c_str(),"w");
} 

void rescaleVelocities(Particle* cube)
{
    /*
     *unsigned int n,i;
     *double vCM[3] = {0,0,0};
     *double l[3] = {0,0,0};
     *double rCM[3]={0,0,0};
     *double omega[3]={0,0,0}; // for calculating angular velocity
     *double dx[3]; // for calculating distance to rCM
     *double invinertia[3][3];
     */

    double vSum = 0;
    for(int i=0;i<N;i++)
        for(int k=0;k<3;k++)
            vSum += cube[i].v[k] * cube[i].v[k];
    double lambda = sqrt(3*(N-1)*Temp/vSum);
    for(int i=0;i<N;i++)
        for(int k=0;k<3;k++)
            cube[i].v[k] *= lambda;
/*
 *    for(i=0;i<3;i++)
 *    {
 *        for(unsigned int j=0;j<3;j++)
 *        {
 *            if(i==j)
 *                invinertia[i][j]=(6./(L*L*N));
 *            else
 *                invinertia[i][j]=0;
 *        }
 *    }
 *
 *    for(n=0;n<N;n++)
 *        for(i=0;i<3;i++)
 *            rCM[i]+=cube[n].r[i];
 *    
 *    for(i=0;i<3;i++)
 *        rCM[i] /= N;
 *        
 *
 *    for(n=0;n<N;n++)
 *        for(i=0;i<3;i++)
 *            vCM[i] += cube[n].v[i];
 *
 *
 *    for(i=0;i<3;i++)
 *        vCM[i] /= N;
 *
 *
 *    for(n=0;n<N;n++)
 *        for(i=0;i<3;i++)
 *            cube[n].v[i] -= vCM[i];
 *
 *    for(n=0;n<N;n++)
 *    {
 *        for(i=0;i<3;i++)
 *            dx[i]=rCM[i]-cube[n].r[i];
 *        l[0]+= (dx[1]*cube[n].v[2]-dx[2]*cube[n].v[1]);
 *        l[1]+= (dx[2]*cube[n].v[0]-dx[0]*cube[n].v[2]);
 *        l[2]+= (dx[0]*cube[n].v[1]-dx[1]*cube[n].v[0]);
 *    }
 *
 *    for(i=0;i<3;i++)
 *        for(unsigned int j=0;j<3;j++)
 *            omega[i]+=invinertia[i][j]*l[j];
 *
 *    for(n=0;n<N;n++)
 *    {
 *        for(i=0;i<3;i++)
 *            dx[i] = rCM[i] - cube[n].r[i];
 *     cube[n].v[0] -= (omega[1]*dx[2]-omega[2]*dx[1]);
 *     cube[n].v[1] -= (omega[2]*dx[0]-omega[0]*dx[2]);
 *     cube[n].v[2] -= (omega[0]*dx[1]-omega[1]*dx[0]);
 *    }
 */
}

void setValues(double temp, double dq, double Eps, double Pressure, double ambienttemp)
{
    Temp = temp;
    dQ = dq;
    eps = Eps;
    P = Pressure;
    AmbientTemp = ambienttemp;
}

/*
 *void calcPressure(Particle* cube,FILE* output)
 *{
 *    double pressure = 0;
 *    double w = 0;
 *    double dist[3];
 *    double distSqd = 0;
 *    
 *    for(int i=0;i<N-1;i++)
 *    {
 *        for(int j=i+1;j<N;j++)
 *        {
 *            for(int k=0;k<3;k++)
 *                dist[k] = cube[i].r[k] - cube[j].r[k];
 *            distSqd = dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2];
 *            w += dist[k]
 *        }
 *    }
 *}
 */

void writeHistory(std::list<Particle*>& gas, int run)
{
    std::list<Particle*>::iterator temp;
    double kinEner;
    double potEner;
    double totEner;
    std::string prefix = "gas/";
    std::string filename;
/*
 *    for(temp = gasHistory.begin();temp != gasHistory.end();temp++)
 *    {
 *        kinEner = 0;
 *        totEner = 0;
 *
 *        for(int i=0;i<3;i++)
 *            kinEner += (*temp)->v[i] * (*temp)->v[i]; 
 *        potEner = (*temp)->potE;
 *        totEner = kinEner+potEner;
 *        kinEner = (*temp)->m * kinEner * 0.5;
 *        (*temp)->kineticEnergy.push_back(std::make_pair((*temp)->name,kinEner));
 *        (*temp)->potentialEnergy.push_back(std::make_pair((*temp)->name,potEner));
 *        (*temp)->totalEnergy.push_back(std::make_pair((*temp)->name,totEner));
 *     //potentialEnergy;
 *     //totalEnergy;
 *    }
 */
    //std::cout << "Length: " << gas.size() << std::endl;
    for(temp = gas.begin();temp != gas.end();temp++)
    {
        FILE* output;
        filename = prefix;
        //filename += "id_";
        filename += numberToString((*temp)->ID);
        filename += ".dat";
        output = fopen(filename.c_str(),"a");
        kinEner = 0;
        totEner = 0;
        for(int i=0;i<3;i++)
            kinEner += (*temp)->v[i] * (*temp)->v[i]; 
        kinEner = (*temp)->m * kinEner * 0.5;
        potEner = (*temp)->potE;
        totEner = kinEner+potEner;
        fprintf(output,"%s \t %lf \t %lf \t %lf \n",((*temp)->name).c_str(),kinEner,potEner,totEner);
        fclose(output);
     }
}

void calculateGasTemperature(std::list<Particle*> gas,FILE* output)
{
    std::list<Particle*>::iterator temp;
    double inTemp = 0;
    double outTemp = 0;
    int inCount = 0;
    int outCount = 0;

    for(temp = gas.begin();temp != gas.end(); temp++)
    {
        if(strcmp((*temp)->name.c_str(),"GasIn") == 0)
        {
            inCount++;
            inTemp += (*temp)->v[0]*(*temp)->v[0]+(*temp)->v[1]*(*temp)->v[1]+(*temp)->v[2]*(*temp)->v[2];
        }
        else if(strcmp((*temp)->name.c_str(),"GasOut") == 0)
        {
            outCount++;
            outTemp += (*temp)->v[0]*(*temp)->v[0]+(*temp)->v[1]*(*temp)->v[1]+(*temp)->v[2]*(*temp)->v[2];
        }
    }

    inTemp = inTemp/(3*inCount);
    outTemp = outTemp/(3*outCount);
    //std::cout << "\tInTemp: " << inTemp << std::endl;
    //std::cout << "\tOutTemp: " << outTemp << std::endl;
    fprintf(output,"%lf\t%lf\t%lf\t%s\n",inTemp,outTemp,inTemp-outTemp,inTemp>outTemp ? "InTemp" : "OutTemp");
}

double calculateEnergies(Particle* cube, std::list<Particle*> gas)
{
    double* energies = new double[3];
    double ePot = 0;
    double eKin = 0;
    double eTot = 0;
    double rij[3] = {0,0,0};
    double rSqd = 0;
    std::list<Particle*>::iterator iter;

    for(int i=0;i<N-1;i++)
    {
        eKin += (cube[i].v[0]*cube[i].v[0]+cube[i].v[1]*cube[i].v[1]+cube[i].v[2]*cube[i].v[2])*0.5;
        for(int j=i+1;j<N;j++)
        {
            for(int m=0;m<3;m++)
                rij[m] = cube[i].r[m] - cube[j].r[m];
            rSqd = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
            if(rSqd < 2.5*2.5)
                ePot += 4*(pow(rSqd,-6)-pow(rSqd,-3));
        }

    }
    for(iter=gas.begin();iter!=gas.end();iter++)
    {
        eKin += ((*iter)->v[0]*(*iter)->v[0]+(*iter)->v[1]*(*iter)->v[1]+(*iter)->v[2]*(*iter)->v[2])*0.5;
        for(int i=0;i<N;i++)
        {
            for(int m=0;m<3;m++)
                rij[m] = cube[i].r[m] - (*iter)->r[m];
            rSqd = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
            if(rSqd < 2.5*2.5)
                ePot += 4*(pow(rSqd,-6)-pow(rSqd,-3));
        }
    }
    
    return eKin+ePot;
        
    
    /*
     *for(int i=0;i<N-1;i++)
     *{
     *    eKin += (cube[i].v[0]*cube[i].v[0]+cube[i].v[1]*cube[i].v[1]+cube[i].v[2]*cube[i].v[2])*0.5;
     *    for(int j=i+1;j<N;j++)
     *    {
     *        for(int m=0;m<3;m++)
     *            rij[m] = cube[i].r[m] - cube[j].r[m];
     *        rSqd = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
     *        if(rSqd < 2.5*2.5)
     *            ePot += 4*(pow(rSqd,-6)-pow(rSqd,-3));
     *    }
     *}
     *eKin += (cube[N].v[0]*cube[N].v[0]+cube[N].v[1]*cube[N].v[1]+cube[N].v[2]*cube[N].v[2])*0.5;
     */

    /*
     *for(iter=gas.begin();iter!=gas.end();iter++)
     *{
     *    eKin += ((*iter)->v[0]*(*iter)->v[0]+(*iter)->v[1]*(*iter)->v[1]+(*iter)->v[2]*(*iter)->v[2])*0.5;
     *    for(int i=0;i<N;i++)
     *    {
     *        for(int m=0;m<3;m++)
     *            rij[m] = cube[i].r[m] - (*iter)->r[m];
     *        rSqd = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
     *        if(rSqd < 2.5*2.5)
     *            ePot += 4*(pow(rSqd,-6)-pow(rSqd,-3));
     *    }
     *}
     */
    /*
     *eTot = eKin+ePot;
     *energies[0] = eKin;
     *energies[1] = ePot;
     *energies[2] = eTot;
     */
    
    //return eTot;
}
double calculateEnergiesTest(Particle* cube, std::list<Particle*> gas)
{
    double* energies = new double[3];
    double ePot = 0;
    double eKin = 0;
    double eTot = 0;
    double rij[3] = {0,0,0};
    double rSqd = 0;
    std::list<Particle*>::iterator iter;

    for(int i=0;i<N-1;i++)
    {
        eKin += (cube[i].v[0]*cube[i].v[0]+cube[i].v[1]*cube[i].v[1]+cube[i].v[2]*cube[i].v[2])*0.5;
        for(int j=i+1;j<N;j++)
        {
            for(int m=0;m<3;m++)
                rij[m] = cube[i].r[m] - cube[j].r[m];
            rSqd = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
            if(rSqd < 2.5*2.5)
                ePot += 4*(pow(rSqd,-6)-pow(rSqd,-3));
        }

    }
    for(iter=gas.begin();iter!=gas.end();iter++)
    {
        eKin += ((*iter)->v[0]*(*iter)->v[0]+(*iter)->v[1]*(*iter)->v[1]+(*iter)->v[2]*(*iter)->v[2])*0.5;
        for(int i=0;i<N;i++)
        {
            for(int m=0;m<3;m++)
                rij[m] = cube[i].r[m] - (*iter)->r[m];
            rSqd = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
            if(rSqd < 2.5*2.5)
                ePot += pow(rSqd,-6);
                //ePot += 4*(pow(rSqd,-6)-pow(rSqd,-3));
        }
    }
    
    return eKin+ePot;
        
    
    /*
     *for(int i=0;i<N-1;i++)
     *{
     *    eKin += (cube[i].v[0]*cube[i].v[0]+cube[i].v[1]*cube[i].v[1]+cube[i].v[2]*cube[i].v[2])*0.5;
     *    for(int j=i+1;j<N;j++)
     *    {
     *        for(int m=0;m<3;m++)
     *            rij[m] = cube[i].r[m] - cube[j].r[m];
     *        rSqd = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
     *        if(rSqd < 2.5*2.5)
     *            ePot += 4*(pow(rSqd,-6)-pow(rSqd,-3));
     *    }
     *}
     *eKin += (cube[N].v[0]*cube[N].v[0]+cube[N].v[1]*cube[N].v[1]+cube[N].v[2]*cube[N].v[2])*0.5;
     */

    /*
     *for(iter=gas.begin();iter!=gas.end();iter++)
     *{
     *    eKin += ((*iter)->v[0]*(*iter)->v[0]+(*iter)->v[1]*(*iter)->v[1]+(*iter)->v[2]*(*iter)->v[2])*0.5;
     *    for(int i=0;i<N;i++)
     *    {
     *        for(int m=0;m<3;m++)
     *            rij[m] = cube[i].r[m] - (*iter)->r[m];
     *        rSqd = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
     *        if(rSqd < 2.5*2.5)
     *            ePot += 4*(pow(rSqd,-6)-pow(rSqd,-3));
     *    }
     *}
     */
    /*
     *eTot = eKin+ePot;
     *energies[0] = eKin;
     *energies[1] = ePot;
     *energies[2] = eTot;
     */
    
    //return eTot;
}

double calculateEnergies(Particle* cube)
{
    double* energies = new double[3];
    double ePot = 0;
    double eKin = 0;
    double eTot = 0;
    double rij[3] = {0,0,0};
    double rSqd = 0;

    for(int i=0;i<N-1;i++)
    {
        eKin += (cube[i].v[0]*cube[i].v[0]+cube[i].v[1]*cube[i].v[1]+cube[i].v[2]*cube[i].v[2])*0.5;
        for(int j=i+1;j<N;j++)
        {
            for(int m=0;m<3;m++)
                rij[m] = cube[i].r[m] - cube[j].r[m];
            rSqd = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
            if(rSqd < 2.5*2.5)
                ePot += 4*(pow(rSqd,-6)-pow(rSqd,-3));
        }

    }

    return eKin+ePot;
        
    
}

void gasStatus(std::list<Particle*> gas)
{
    unsigned int gIn = 0;
    unsigned int gChng = 0;
    unsigned int gOut = 0;
    std::list<Particle*>::iterator gasIter;
    for(gasIter = gas.begin();gasIter!=gas.end();gasIter++)
    {
        if(strcmp((*gasIter)->name.c_str(),"GasIn") == 0)
            gIn++;
        if(strcmp((*gasIter)->name.c_str(),"GasChange") == 0)
            gChng++;
        if(strcmp((*gasIter)->name.c_str(),"GasOut") == 0)
            gOut++;
    }

    std::cout << "In: " << gIn << "\t" << "Change: " << gChng << "\t" << "Out: " << gOut << std::endl;

}

void verletBaroAccelerations(Particle* cube, std::list<Particle*> gas)
{
    double fVerlet = 0;
    double fBaro = 0;
    double fHarm[3] = {0,0,0};
    double rij[3] = {0,0,0};
    double rSqd = 0;
    double rCOM[3] = {0,0,0};
    double dist[3] = {0,0,0};

    // INITIALIZE ACCELERATIONS

    std::list<Particle*>::iterator gasIter;
    for(unsigned int i=0;i<N;i++)
        for(unsigned int m=0;m<3;m++)
        {
            cube[i].aOld[m] = cube[i].a[m];
            cube[i].a[m] = 0.;        
        }
    for(gasIter=gas.begin();gasIter != gas.end();gasIter++)
       for(unsigned int m=0;m<3;m++)
          (*gasIter)->a[m] = 0.; 

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
        for(gasIter=gas.begin();gasIter != gas.end();gasIter++)
        {
            rSqd = 0;
            for(int m=0;m<3;m++)
              rij[m] = (*gasIter)->r[m] - cube[i].r[m];
            rSqd = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
            if(rSqd <= rCut2)
            {
                fBaro = 12. * pow(rSqd,-6.);
                for(int m=0;m<3;m++)
                {
                    (*gasIter)->a[m] += rij[m] * fBaro;
                    cube[i].a[m] -= rij[m] * fBaro;
                }
            }
        }
    }
    // HARMONIC TRAP
    for(unsigned int i=0;i<N;i++)
        for(unsigned int k=0;k<3;k++)
            rCOM[k] += cube[i].r[k]; 

    for(unsigned int k=0;k<3;k++)
        rCOM[k] /= N;

    for(unsigned int k=0;k<3;k++)
       dist[k] = center[k]-rCOM[k]; 

    //for(unsigned int k=0;k<3;k++)
       //fHarm[k] += k*dist[k]; 

    for(unsigned int i=0;i<N;i++)
        for(unsigned int k=0;k<3;k++)
            cube[i].a[k] += k*dist[k];
    

}
void verletBaro(Particle* cube, std::list<Particle*>& gas)
{
    InitBarostatFull(gas);
    std::list<Particle*>::iterator gasIter;
    for(unsigned int i=0;i<N;i++)
    {
        for(unsigned int k=0;k<3;k++)
        {
            cube[i].r[k] +=cube[i].v[k]*dt + 0.5 * cube[i].a[k] * dt * dt;
            cube[i].v[k] += 0.5 * cube[i].a[k] * dt;
        }         
    }
    for(gasIter=gas.begin();gasIter != gas.end();gasIter++)
    {
        for(unsigned int k=0;k<3;k++)
        {
            (*gasIter)->r[k] +=(*gasIter)->v[k]*dt + 0.5 * (*gasIter)->a[k] * dt * dt / (*gasIter)->m;
            (*gasIter)->v[k] += 0.5 * (*gasIter)->a[k] * dt / (*gasIter)->m;
        }         
    }

    CheckBoundariesNew(cube,gas);
    verletBaroAccelerations(cube,gas);

    for(unsigned int i=0;i<N;i++)
        for(unsigned int k=0;k<3;k++)
            cube[i].v[k] += 0.5 * cube[i].a[k] * dt;
    
    for(gasIter=gas.begin();gasIter != gas.end();gasIter++)
        for(unsigned int k=0;k<3;k++)
            (*gasIter)->v[k] += 0.5 * (*gasIter)->a[k] * dt / (*gasIter)->m;
}

void eHEXBaro(Particle* cube, std::list<Particle*>& gas)
{
	double eps[N][3];
	double K=0;
	double eta[N][3];
	double xi=0;
	double sumforces[3]={0,0,0};
	double sumvel=0; 
    unsigned int i,n;

    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            K += (cube[n].v[i]*cube[n].v[i]);
    K = K/2.;
    xi = sqrt(1+dQ/K);
    
    InitBarostatFull(gas);
    std::list<Particle*>::iterator gasIter;

    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            cube[n].vnew[i] = xi*cube[n].v[i]; //(19a)
    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            cube[n].vhalf[i] = cube[n].vnew[i]+dt*0.5*cube[n].a[i]; //(19b)
    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            cube[n].rnew[i] = cube[n].r[i]+dt*cube[n].vhalf[i]; //(19c)
    for(n=0;n<N;n++)
    {
        for(i=0;i<3;i++)
        {
            cube[n].v[i] = cube[n].vhalf[i];
            cube[n].r[i] = cube[n].rnew[i];
        }
    }
    /*
     *for(unsigned int i=0;i<N;i++)
     *{
     *    for(unsigned int k=0;k<3;k++)
     *    {
     *        cube[i].r[k] +=cube[i].v[k]*dt + 0.5 * cube[i].a[k] * dt * dt;
     *        cube[i].v[k] += 0.5 * cube[i].a[k] * dt;
     *    }         
     *}
     */
    for(gasIter=gas.begin();gasIter != gas.end();gasIter++)
    {
        for(unsigned int k=0;k<3;k++)
        {
            (*gasIter)->r[k] +=(*gasIter)->v[k]*dt + 0.5 * (*gasIter)->a[k] * dt * dt / (*gasIter)->m;
            (*gasIter)->v[k] += 0.5 * (*gasIter)->a[k] * dt / (*gasIter)->m;
        }         
    }

    CheckBoundariesNew(cube,gas);
    verletBaroAccelerations(cube,gas);
    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            cube[n].vnew[i] = cube[n].vhalf[i] + dt*0.5*cube[n].a[i]; //(19e)
    K = 0; 
    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            K += (cube[n].vnew[i]*cube[n].vnew[i]);
    K = K/2.;
    xi = sqrt(1+dQ/K); //\bar{xi} from (19f)
    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            cube[n].v[i] = xi*cube[n].vnew[i]; //(19f)

	for(n=0;n<N;n++)
		for(i=0;i<3;i++)
			eta[n][i]=(FG/(2.*K))*cube[n].vnew[i];
    sumvel=0;
    for(i=0;i<3;i++)
        sumforces[i]=0;
	for(n=0;n<N;n++)
		for(i=0;i<3;i++)
		{
			sumforces[i]+=cube[n].a[i];
			sumvel+=cube[n].vnew[i]*cube[n].a[i];
		}
	for(n=0;n<N;n++)
		for(i=0;i<3;i++)
            eps[n][i]=(eta[n][i]/K)*((FG/48.)+(1./6.)*sumvel)-(FG/(12.*K))*
            (cube[n].a[i]-(sumforces[i])); //(20)
    for(n=0;n<N;n++)
        for(i=0;i<3;i++)
            cube[n].r[i] = cube[n].rnew[i] - dt*dt*dt*eps[n][i]; //(19g)

    /*
     *for(unsigned int i=0;i<N;i++)
     *    for(unsigned int k=0;k<3;k++)
     *        cube[i].v[k] += 0.5 * cube[i].a[k] * dt;
     */
    
    for(gasIter=gas.begin();gasIter != gas.end();gasIter++)
        for(unsigned int k=0;k<3;k++)
            (*gasIter)->v[k] += 0.5 * (*gasIter)->a[k] * dt / (*gasIter)->m;
}
