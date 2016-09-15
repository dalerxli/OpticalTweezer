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


void addToList(std::list<Particle*>& particles,std::list<Particle*>& particlesHistory)
{
    int i = rand() % 100;
    std::cout << i << std::endl;
    if(particles.size() != 0)
    {
        std::list<Particle*>::iterator temp = particles.end();
        for(int k=0;k<i;k++)
            particles.push_back(new Particle());
        for(temp;temp != particles.end();temp++)
        {
            particlesHistory.push_back(*temp);
        }
    }
    else
    {
        std::list<Particle*>::iterator temp;
        for(int k=0;k<i;k++)
            particles.push_back(new Particle());
        for(temp=particles.begin();temp != particles.end();temp++)
        {
            particlesHistory.push_back(*temp);
        }
    }


    
}


int main(int argc, char** argv)
{
    std::list<Particle*> particles;
    std::list<Particle*> particlesHistory;
    std::list<Particle*>::iterator pIter;
    std::list<Particle*>::iterator pHistIter;
    for(int i=0;i<5;i++)
        addToList(particles,particlesHistory);
    pHistIter = particlesHistory.begin();
    for(pIter=particles.begin();pIter != particles.end();pIter++)
    {
        std::cout << (*pIter)->ID << "\t";
        std::cout << (*pHistIter)->ID << std::endl;
        pHistIter++;
    }
    
	return 0;
}
