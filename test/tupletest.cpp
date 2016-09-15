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

int main(int argc, char** argv)
{
    std::vector<std::pair<int,double> > test1;
    test1.push_back(std::make_pair(1,3.14));
    std::cout << test1[0].second<< std::endl;

    std::pair<std::string,double> test2; 
    test2 = std::make_pair("test!",3.14);
    std::cout << "(" << test2.first << "," << test2.second << ")" << std::endl;


	return 0;

    
}
