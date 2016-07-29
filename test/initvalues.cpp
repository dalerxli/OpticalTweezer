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
    std::string dir = "test";
    chdir(dir.c_str());
    FILE* initvalues = fopen("initvalues.dat","wb");
    double T;
    double dQ;
    std::cout << "\tT = ";
    std::cin >> T;  
    std::cout << "\tdQ = ";
    std::cin >> dQ;  
    fwrite(&T,sizeof(double),1,initvalues);
    fwrite(&dQ,sizeof(double),1,initvalues);
    
    fclose(initvalues);
    return 0;
}
