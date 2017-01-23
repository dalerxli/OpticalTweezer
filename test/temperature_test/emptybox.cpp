#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include "../../include/functions.hpp"
#include "../../include/globals.hpp"
#include "../../include/classes.hpp"

int main(int argc, char** argv)
{
    setValues(0.2,0.07,0.2,4,0.1);
    Particle *cube = new Particle[N];
    InitPositions(cube);
    InitVelocities(cube);
    ComputeAccelerations(cube);

    FILE* output = fopen("emptybox_gas_data.dat","w");
    int run;

    std::list<Particle*> gas;
    std::list<Particle*>::iterator temp;

    ComputeSoftSphereNoCube(gas);
    for(run=0;run<50000;run++)
    {
        if(run%500==0)
            printf("(MEASURE) Zeitschritt %d\n",run);
        BarostatNoCube(gas);
        std::cout << gas.size() << std::endl;
/*
 *        double inTemp = 0;
 *        double outTemp = 0;
 *        int inCount = 0;
 *        int outCount = 0;
 *
 *        for(temp = gas.begin();temp != gas.end(); temp++)
 *        {
 *            if(strcmp((*temp)->name.c_str(),"GasIn") == 0)
 *            {
 *                inCount++;
 *                inTemp += (*temp)->v[0]*(*temp)->v[0]+(*temp)->v[1]*(*temp)->v[1]+(*temp)->v[2]*(*temp)->v[2];
 *            }
 *            else if(strcmp((*temp)->name.c_str(),"GasOut") == 0)
 *            {
 *                outCount++;
 *                outTemp += (*temp)->v[0]*(*temp)->v[0]+(*temp)->v[1]*(*temp)->v[1]+(*temp)->v[2]*(*temp)->v[2];
 *            }
 *        }
 *
 *        inTemp = inTemp/(3*inCount);
 *        outTemp = outTemp/(3*outCount);
 *        fprintf(output,"%lf\t%lf\t%lf\t%s\n",inTemp,outTemp,inTemp-outTemp,inTemp>outTemp ? "InTemp" : "OutTemp");
 */
    }

    fclose(output);

	return 0;
}
