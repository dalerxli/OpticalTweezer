#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <vector>

int main(int argc, char** argv)
{
    std::string test1 = "GasGone";
    std::string test2 = "GbtGone";
    std::string test3 = "GasGone";
    std::string test4 = "GasGone";
    std::string test5 = "GasGone";
    std::string test6 = "GasGone";
    std::string test7 = "GasGone";

    int a1 = strcmp(test1.c_str(),test2.c_str());    
    int a2 = strcmp(test1.c_str(),test3.c_str());    
    int a3 = strcmp(test1.c_str(),test4.c_str());    
    int a4 = strcmp(test1.c_str(),test5.c_str());    
    int a5 = strcmp(test1.c_str(),test6.c_str());    
    int a6 = strcmp(test1.c_str(),test7.c_str());    

    std::cout << a1 << std::endl;
    std::cout << a2 << std::endl;
    std::cout << a3 << std::endl;
    std::cout << a4 << std::endl;
    std::cout << a5 << std::endl;
    std::cout << a6 << std::endl;

    return 0;
}
