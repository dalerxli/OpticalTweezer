#include "../include/globals.hpp"
#include "../include/classes.hpp"
#include "../include/functions.hpp"

int main()
{
    Particle* particles = new Particle[N];
    std::ifstream inFile ("../output/equidata.dat");
    
    for(int i=0;i<N;i++)
    {
        inFile >> particles[i].r[0];
        inFile >> particles[i].r[1];
        inFile >> particles[i].r[2];
        inFile >> particles[i].v[0];
        inFile >> particles[i].v[1];
        inFile >> particles[i].v[2];
        inFile >> particles[i].a[0];
        inFile >> particles[i].a[1];
        inFile >> particles[i].a[2];
    }
    for(int i=0;i<N;i++)
    {
        std::cout << std::setprecision(9)<< particles[i].r[0] << "\t";
        std::cout << std::setprecision(9)<< particles[i].r[1] << "\t";
        std::cout << std::setprecision(9)<< particles[i].r[2] << "\t";
        std::cout << std::setprecision(9)<< particles[i].v[0] << "\t";
        std::cout << std::setprecision(9)<< particles[i].v[1] << "\t";
        std::cout << std::setprecision(9)<< particles[i].v[2] << "\t";
        std::cout << std::setprecision(9)<< particles[i].a[0] << "\t";
        std::cout << std::setprecision(9)<< particles[i].a[1] << "\t";
        std::cout << std::setprecision(9)<< particles[i].a[2] << "\t";
        std::cout << std::endl;
    }
    inFile.close();
    return 0;
}
