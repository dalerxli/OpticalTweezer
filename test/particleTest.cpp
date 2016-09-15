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
    Particle* p1 = new Particle();
    Particle  *p2,*p3;
    std::cout << p1->ID << std::endl;
    for(int i=0;i<10;i++)
    {
        p2 = new Particle();
        std::cout << p2->ID << std::endl;
        delete p2;
    }
    p2 = new Particle();
    p3 = new Particle();

    delete p1;
    delete p2;
    delete p3;

    /*
     *for(int i=0;i<N;i++)
     *{
     *    for(int j=0;j<N;j++)
     *    {
     *        delete [] Forces[i][j];
     *        delete [] Distances[i][j];
     *    }
     *}
     */
    for(int i=0;i<N;i++)
    {
        delete [] Forces[i];
        delete [] Distances[i];
    }
    delete [] Forces;
    delete [] Distances;
//double*** Forces = new double**[N];
//double*** Distances = new double**[N];
	return 0;
}
