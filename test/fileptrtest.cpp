#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include "../include/globals.hpp"
#include "../include/classes.hpp"
#include "../include/functions.hpp"

void writeFile();

int main(int argc, char** argv)
{
    for(int i=0;i<350;i++)
    {
        std::cout << i << std::endl;
        writeFile();
    }
	return 0;
}

void writeFile()
{
    std::string prefix = "fileptrtest/";
    std::string filename;
    for(int i=0;i<10;i++)
    {
        FILE* testout;
        filename = prefix;
        filename += numberToString(i);
        filename += ".dat";
        testout = fopen(filename.c_str(),"a");
        fprintf(testout,"Hello world from %d\n",i);
        fclose(testout);
    }
    //free(testout);
    //delete testout;
}
