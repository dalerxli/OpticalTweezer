#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include "../include/globals.hpp"
#include "../include/classes.hpp"
#include "../include/functions.hpp"

int main(int argc, char** argv)
{
    FILE* initvalues = fopen("initvalues.dat","rb");
    double T;
    double dQ;
    fread(&T,sizeof(double),1,initvalues);
    fread(&dQ,sizeof(double),1,initvalues);

    std::cout << T << std::endl;
    std::cout << dQ << std::endl;
    
    fclose(initvalues);
    return 0;
}
