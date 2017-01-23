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
    Particle *particles = new Particle[N];
    InitPositions(particles);
    double rcom[3] = {0,0,0};
    for(int i=0;i<N;i++)
        for(int k=0;k<3;k++)
            rcom[k] += particles[i].r[k];
    for(int k=0;k<3;k++)
        rcom[k] /= N;

    std::cout << "rCOM: " << rcom[0] << ", " << rcom[1] << ", " << rcom[2] << std::endl;
    std::cout << "L/2: " << L/2. << ", " << L/2. << ", " << L/2. << std::endl;

	return 0;
}
