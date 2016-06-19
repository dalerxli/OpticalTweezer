#include "../include/globals.hpp"
#include "../include/classes.hpp"
#include "../include/functions.hpp"

int main(int argc,char** argv)
{
    system("rm output/combined/*");
    int a;
    int run=0;
    gsl_rng_set(r,98);
    FILE* output = fopen("output/box_test.xyz","w");
    FILE* pressure = fopen("output/pressure.dat","w");
    FILE* baroOut = fopen("output/barostat.dat","w");
    FILE* combinedOut = fopen("output/combined.xyz","w");
    FILE* particledata = fopen("output/particledata.dat","w");
    FILE* tempout = fopen("output/temperature.dat","w");
    FILE* cmData = fopen("output/cmdata.dat","w");
    FILE* equiData = fopen("output/equidata.dat","w");
    FILE* readData = fopen("output/equidata.dat","r");
    FILE* totErgData = fopen("output/energy.dat","w");
    FILE* equidata = fopen("output/readwritetest.dat","w");
    FILE* equidata2 = fopen("output/readwritetest2.dat","w");
    FILE* surfacedata = fopen("output/surfacedata.xyz","w");
    FILE* particleTracker = fopen("output/combined/id_1400.xyz","w");
    fclose(combinedOut);
    Particle *cube = new Particle[N];
    std::list<Particle*> gas;
    double* rCM = new double[3];
    double* vCM = new double[3];
    double energy = 0;
    InitPositions(cube);
    fclose(surfacedata);
    calcCM(cube,rCMStart,vCM);
    calcCM(cube,rCM,vCM);
    fprintf(output,"%d\nCube\n",N);
    InitVelocities(cube);
    fprintf(pressure,"%lf\n",Pressure(cube));
    ComputeAccelerations(cube);
/*
 *    for(run=0;run<10000;run++)
 *    {
 *            if(run%200==0)
 *                printf("Zeitschritt %d\n",run);
 *            VelocityVerlet(cube,0,output);
 *            if(run%100 == 0)
 *            {
 *                calcTemp(cube,tempout); 
 *                //energy = totErg(cube);
 *                //fprintf(totErgData,"%d\t%lf\n",run,energy);
 *            }
 *    }
 *
 *    writePositions(cube,"output/states/LJequilibriumSur.dat");
 */

    /*
     *for(run;run<50000;run++)
     *{
     *        if(run%200==0)
     *            printf("Zeitschritt %d\n",run);
     *        fprintf(output,"864\nNew Frame\n");
     *        eHEX(cube,output);
     *        fprintf(pressure,"%lf\n",Pressure(cube));
     *        calcTemp(cube,tempout); 
     *        calcCM(cube,rCM,vCM);
     *        fprintf(cmData,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
     *                rCM[0],rCM[1],rCM[2],sqrt(rCM[0]*rCM[0]+rCM[1]*rCM[1]+rCM[2]*rCM[2]),vCM[0],vCM[1],vCM[2],sqrt(vCM[0]*vCM[0]+vCM[1]*vCM[1]+vCM[2]*vCM[2]));
     *}
     */
    
    readPositions(cube,"output/states/LJequilibriumSur.dat");
    //readPositions(cube,"output/states/eHEXEquilibrium.dat");
        
    ComputeSoftSphere(gas,cube);

    for(run = 0;run<20000;run++)
    {
        if(run%100==0)
            printf("Zeitschritt %d\n",run);
        eHEX(cube,output);
        Barostat(cube,gas);
        harmonicTrap(rCM,vCM,rCMStart,cube);
        if(run%10==0)
        {
            trackParticle(cube,gas,1400,particleTracker);
            GenerateOutput(cube,gas,run);
        }
        if(run%10==0)
            calcTemp(cube,tempout); 
    }

    //writePositions(cube,"output/states/eHEXEquilibrium.dat");




    fclose(pressure);
    fclose(output);
    fclose(baroOut);
    fclose(particledata);
    fclose(tempout);
    fclose(cmData);
    fclose(equiData);
    fclose(totErgData);
    delete [] cube;
    return 0;
}
