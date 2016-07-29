#include "../include/globals.hpp"
#include "../include/classes.hpp"
#include "../include/functions.hpp"

int main(int argc,char** argv)
{
    /*
     *if(argc > 1)
     *    input = argv[1];
     *else
     *    strcpy(input,"");
     *std::string makedir = "mkdir output/combined_";
     *std::string dirname(input);
     *makedir += dirname;
     *system(makedir.c_str());
     *system("rm output/combined/*");
     */
    std::string folderName = "output/runs/";
    folderName += DateToString();
    int dir = mkdir(folderName.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(dir == -1)
    {
        folderName += "_1";
        mkdir(folderName.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    chdir(folderName.c_str()); 
    std::ofstream params;
    params.open("parameters.txt");
    params << rho;
    params << Temp;
    params << dt;
    params << eps;
    params << dQ;
    params.close();
/*
 *double FG;
 *double* rCMStart = new double[3];
 *double rho = 1.1;
 *double L = pow(N/rho,1.0/3);
 *double rCutOff = 2.5;
 *double dt = 0.01;
 *double eps = 0.5;
 *double dQ=0.050;
 *char* input;
        double mu;
        //double P=0.0024;
        //This is the calculated atmospheric pressure with LJ units
        double P=4;
        double k=1.0;
        double T=0.9;
 */
    int a;
    int run=0;
    gsl_rng_set(r,98);
    /*
     *FILE* output = fopen("output/box_test.xyz","w");
     *FILE* pressure = fopen("output/pressure.dat","w");
     *FILE* baroOut = fopen("output/barostat.dat","w");
     *FILE* combinedOut = fopen("output/combined.xyz","w");
     *FILE* particledata = fopen("output/particledata.dat","w");
     *FILE* tempout = fopen("output/temperature.dat","w");
     *FILE* cmData = fopen("output/cmdata.dat","w");
     *FILE* equiData = fopen("output/equidata.dat","w");
     *FILE* readData = fopen("output/equidata.dat","r");
     *FILE* totErgData = fopen("output/energy.dat","w");
     *FILE* equidata = fopen("output/readwritetest.dat","w");
     *FILE* equidata2 = fopen("output/readwritetest2.dat","w");
     *FILE* surfacedata = fopen("output/surfacedata.xyz","w");
     *FILE* particleTracker = fopen("output/combined/id_1400.xyz","w");
     */
    //FILE* output = fopen("box_test.xyz","w");
    //FILE* pressure = fopen("output/pressure.dat","w");
    //FILE* baroOut = fopen("output/barostat.dat","w");
    //FILE* combinedOut = fopen("output/combined.xyz","w");
    //FILE* particledata = fopen("output/particledata.dat","w");
    FILE* tempout = fopen("temperature_internal.dat","w");
    FILE* comData = fopen("comdata.dat","w");
    //FILE* cmData = fopen("output/cmdata.dat","w");
    //FILE* equiData = fopen("output/equidata.dat","w");
    //FILE* readData = fopen("output/equidata.dat","r");
    //FILE* totErgData = fopen("output/energy.dat","w");
    //FILE* equidata = fopen("output/readwritetest.dat","w");
    //FILE* equidata2 = fopen("output/readwritetest2.dat","w");
    //FILE* surfacedata = fopen("output/surfacedata.xyz","w");
    //FILE* particleTracker = fopen("output/combined/id_1400.xyz","w");
    //fclose(combinedOut);
    std::cout << "works!" << std::endl;
    Particle *cube = new Particle[N];
    std::list<Particle*> gas;
    double* rCM = new double[3];
    double* rCMtemp = new double[3];
    double* vCM = new double[3];
    double energy = 0;
    InitPositions(cube);
    //fclose(surfacedata);
    calcCM(cube,rCMStart,vCM);
    calcCM(cube,rCM,vCM);
    InitVelocities(cube);
    //fprintf(pressure,"%lf\n",Pressure(cube));
    ComputeAccelerations(cube);
    for(run=0;run<10000;run++)
    {
            if(run%200==0)
                printf("Zeitschritt %d\n",run);
            VelocityVerlet(cube,0,NULL);
            if(run%100 == 0)
                calcTemp(cube,tempout); 
    }

    //writePositions(cube,"output/states/LJequilibriumSur.dat");

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
    
    //readPositions(cube,"output/states/LJequilibriumSur.dat");
    //readPositions(cube,"output/states/eHEXEquilibrium.dat");
        
    ComputeSoftSphere(gas,cube);

    for(run = 0;run<20000;run++)
    {
        if(run%100==0)
            printf("Zeitschritt %d\n",run);
        eHEX(cube);
        Barostat(cube,gas);
        //harmonicTrap(rCM,vCM,rCMStart,cube);
        if(run%10==0)
        {
            //trackParticle(cube,gas,1400,particleTracker);
            GenerateOutput(cube,gas,run);
            calcCM(cube,rCMtemp,comData);
        }
        if(run%10==0)
            calcTemp(cube,tempout); 
    }

    //writePositions(cube,"output/states/eHEXEquilibrium.dat");




    /*
     *fclose(pressure);
     *fclose(output);
     *fclose(baroOut);
     *fclose(particledata);
     *fclose(tempout);
     *fclose(cmData);
     *fclose(equiData);
     *fclose(totErgData);
     */
    delete [] cube;
    return 0;
}
