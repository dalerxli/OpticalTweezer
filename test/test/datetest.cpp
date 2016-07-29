#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include "../../include/globals.hpp"
#include "../../include/functions.hpp"
#include "../../include/classes.hpp"
#include <sys/stat.h>

int main(int argc, char** argv)
{
    std::string test = DateToString();
    std::cout << test << std::endl;
    int dir = mkdir(test.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    chdir(test.c_str()); 
    system("pwd");
    FILE* testout1 = fopen("testfile1.txt","w");
    fprintf(testout1,"Test!\n");
    fclose(testout1);
    chdir("../");
    sleep(65);
    //std::cout << dir << std::endl;
    test = DateToString();
    dir = mkdir(test.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(dir == -1)
    {
        test += "_1";
        mkdir(test.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    std::cout << test << std::endl;
    chdir(test.c_str()); 
    system("pwd");
    FILE* testout2 = fopen("testfile2.txt","w");
    fprintf(testout2,"Test!\n");
    fclose(testout2);
    chdir("../");

    return 0;
}
