#include "../include/globals.hpp"
#include "../include/classes.hpp"
#include "../include/functions.hpp"

void mainLoop();


int main(int argc,char** argv)
{
    //int returnValue;
/*
 *void setValues(double temp, double dq, double Eps, double Pressure, double ambienttemp)
 *double rho = 1.1;
 *double Temp = 0.2;
 *double AmbientTemp = 0.9;
 *double P= 4.0;
 *double L = pow(N/rho,1.0/3);
 *double rCutOff = 2.5;
 *double dt = 0.01;
 *double eps = 0.5;
 *double dQ=0.050;
 */
    for(unsigned int i=0;i<N;i++)
    {
        Forces[i] = new double*[N];
        Distances[i] = new double*[N];
    }
    for(unsigned int i=0;i<N;i++)
        for(unsigned int j=0;j<N;j++)
        {
            Forces[i][j] = new double[3];
            Distances[i][j] = new double[3];
        }

    double t = 0.2;
    double ambient = 0.06;
    double press = 0.5;
    double q = 0.04;
    while(q < 0.1)
    {
        ambient = 0.06;
        while(ambient < 0.1)
        {
            setValues(t,q,0.2,press,ambient);
            std::cout << "==================================================" << std::endl;
            std::cout << "t: " << t << std::endl;
            std::cout << "p: " << press << std::endl;
            std::cout << "q: " << q << std::endl;
            std::cout << "at: " << ambient << std::endl;
            std::cout << "==================================================" << std::endl;
            mainLoop();
            ambient += 0.01;
        }
        q += 0.01;
    }
    //setValues(0.2,0.04,0.2,0.8,0.06);
    //mainLoop();
    return 0;
}


void mainLoop() {
    /*
     *if(argc > 1)
     *    input = argv[1];
     *else
     *    strcpy(input,"");
     *std::string makedir = "mkdir output/combined_";
     *std::string dirname(input);
     *makedir += dirname;
     *system(makedir.c_str());
     *system("rm output/combined/ *");
     */
    std::string folderName = "output/runs/";
    folderName += DateToString();
    int dir = mkdir(folderName.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    int i=1;
    while(dir == -1)
    {
        std::string newFolderName = folderName + "_" + numberToString(i);
        //folderName += "_";
        //folderName += numberToString(i);
        dir = mkdir(newFolderName.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if(dir == 0)
            folderName = newFolderName;
        i++;
    }
    chdir(folderName.c_str()); 
    system("mkdir combined");
    system("mkdir gas");
    std::ofstream params;
    params.open("parameters.txt");
    params << "rho = " << rho << std::endl;
    params << "Temp = " << Temp << std::endl;
    params << "dt = " << dt << std::endl;
    params << "eps = " << eps << std::endl;
    params << "P = " << P << std::endl;
    params << "dQ = " << dQ << std::endl;
    params << "AmbientTemp = " << AmbientTemp << std::endl;
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
    gsl_histogram_set_ranges_uniform(gas_in,-0.1,3);
    gsl_histogram_set_ranges_uniform(gas_out,-0.1,3);
    gsl_histogram_set_ranges_uniform(gas_real_in,-0.1,3);
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
    FILE* tempout = fopen("temperature_internal.dat","w");
    FILE* COMtempout = fopen("temperature_com.dat","w");
    FILE* comData = fopen("comdata.dat","w");
    FILE* vCOMData = fopen("vCOMData.dat","w");
    FILE* gasInData = fopen("histogram_in.dat","w");
    FILE* gasOutData = fopen("histogram_out.dat","w");
    FILE* gasRealInData = fopen("histogram_real_in.dat","w");
    FILE* gasTempData = fopen("gasTempData.dat","w");
    Particle *cube = new Particle[N];
    std::list<Particle*> gas;
    std::list<Particle*> gasHistory;
    double* rCM = new double[3];
    double* rCMtemp = new double[3];
    double* vCM = new double[3];
    double energy = 0;
    InitPositions(cube);
    calcCM(cube,rCMStart,vCM);
    calcCM(cube,rCM,vCM);
    InitVelocities(cube);
    ComputeAccelerations(cube);
    for(run=0;run<10000;run++)
    {
            if(run%500==0)
                printf("(INIT) Zeitschritt %d\n",run);
            VelocityVerlet(cube,0,NULL);
            if(run%100 == 0)
            {
                rescaleVelocities(cube);
                calcTemp(cube,tempout); 
                //calcCM(cube,rCMtemp,comData);
                for(int k=0;k<3;k++)
                {
                    rCM[k] = 0;
                    vCM[k] = 0;
                }
                for(unsigned int j=0;j<N;j++)
                    for(int k=0;k<3;k++)
                    {
                        rCM[k] += cube[j].r[k];
                        vCM[k] += cube[j].v[k];
                    }
                for(int k=0;k<3;k++)
                {
                    rCM[k] = rCM[k]/(N*1.0);
                    vCM[k] = vCM[k]/(N*1.0);
                }
                fprintf(comData,"%lf\t%lf\t%lf\n",rCM[0],rCM[1],rCM[2]);
                fprintf(vCOMData,"%lf\t%lf\t%lf\n",vCM[0],vCM[1],vCM[2]);
            }
            //if(run%400 == 0)
                //GenerateOutput(cube,gas,run);
    }

    //writePositions(cube,"../../states/LJequilibriumSur.dat");

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
    
    //readPositions(cube,"../../states/LJequilibriumSur.dat");
    //readPositions(cube,"output/states/eHEXEquilibrium.dat");
        
    ComputeSoftSphere(gas,cube);

    for(run = 0;run<80000;run++)
    {
        if(run%500==0)
        {
            printf("(MEASURE) Zeitschritt %d - Number of Gas particles: %d\n",run,gas.size());
        }
        eHEX(cube);
        BarostatNew(cube,gas);
        //BarostatNew(cube,gas,gasHistory);
        //std::cout << "works!" << std::endl;
        harmonicTrap(rCM,vCM,rCMStart,cube);
        //std::cout << "works!" << std::endl;
        if(run%400==0)
        {
            //trackParticle(cube,gas,1400,particleTracker);
            GenerateOutput(cube,gas,run+10000);
            //std::cout << "works!" << std::endl;
            calcCM(cube,rCMtemp,comData);
        }
        if(run%100==0)
        {
            calcTemp(cube,tempout); 
            //calcCOMTemp(vCM,COMtempout);
            //fprintf(vCOMData,"%lf\t%lf\t%lf\n",vCM[0],vCM[1],vCM[2]);
            calculateGasTemperature(gas,gasTempData);
            for(int k=0;k<3;k++)
            {
                rCM[k] = 0;
                vCM[k] = 0;
            }
            for(unsigned int j=0;j<N;j++)
                for(int k=0;k<3;k++)
                {
                    rCM[k] += cube[j].r[k];
                    vCM[k] += cube[j].v[k];
                }
            for(int k=0;k<3;k++)
            {
                rCM[k] = rCM[k]/(N*1.0);
                vCM[k] = vCM[k]/(N*1.0);
            }
            fprintf(comData,"%lf\t%lf\t%lf\n",rCM[0],rCM[1],rCM[2]);
            fprintf(vCOMData,"%lf\t%lf\t%lf\n",vCM[0],vCM[1],vCM[2]);
        }
        //std::cout << "works!" << std::endl;
        //writeHistory(gas,run);
    }
    /*
     *FILE* pressure = fopen("inst_pressure.dat","w");
     *for(unsigned int i=0;i<virial.size();i++)
     *    fprintf(pressure,"%lf\n",virial[i]);
     *fclose(pressure);
     */

    gsl_histogram_fprintf(gasInData,gas_in,"%f","%f");
    gsl_histogram_fprintf(gasOutData,gas_out,"%f","%f");
    gsl_histogram_fprintf(gasRealInData,gas_real_in,"%f","%f");
/*
 *    FILE* id_test = fopen("id_test.dat","w");
 *    std::list<Particle*>::iterator p_iter = gasHistory.begin();
 *    for(int i=0;i<10;i++)
 *        ++p_iter;
 *    for(int i=0;i<(*p_iter)->kineticEnergy.size();i++)
 *        fprintf(id_test,"%s \t %lf \n",((*p_iter)->kineticEnergy[i].first).c_str(),(*p_iter)->kineticEnergy[i].second);
 *     //std::vector<std::pair<std::string,double> >  kineticEnergy;
 *
 *    fclose(id_test);
 */
    chdir("../../../");
    delete [] cube;
    gas.clear(); 
    g_ID = 1;
	EHEX_FLAG = false;
    fclose(tempout);
    fclose(comData);
    fclose(COMtempout);
    fclose(vCOMData);
    fclose(gasInData);
    fclose(gasOutData);
    fclose(gasRealInData);
    fclose(gasTempData);
    for(int i=0;i<N;i++)
    {
        delete [] Forces[i];
        delete [] Distances[i];
    }
    delete [] Forces;
    delete [] Distances;
}
