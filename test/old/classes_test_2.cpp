#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <list>
#include <iterator> 

int g_ID=1;
const gsl_rng_type *T = gsl_rng_taus;
gsl_rng *r = gsl_rng_alloc(T);

class Particle
{
	public:
	 double x,y,z;
	 double vx,vy,vz;
	 int ID;
	 Particle();
};

Particle::Particle()
{
				x = gsl_rng_uniform(r);
				y = gsl_rng_uniform(r);
				z = gsl_rng_uniform(r);
				vx = gsl_rng_uniform(r);
				vy = gsl_rng_uniform(r);
				vz = gsl_rng_uniform(r);
				ID = g_ID;
				g_ID++;
}

void PrintListElement(std::list<Particle*>::const_iterator particleIterator)
{
				std::cout << (*particleIterator)->ID << ",("<< 
								(*particleIterator)->x << ","<<(*particleIterator)->y <<","<<
								(*particleIterator)->z <<")"<< std::endl;	
}

int main()
{
				std::list<Particle*> particleList;
				//Particle *parList = new Particle[100];
				gsl_rng_set(r,132);

				for(int i=0;i<100;i++)
				{
						particleList.push_back(new Particle());
				}

				std::list<Particle*>::iterator particleIterator;
				particleIterator = particleList.begin();
				PrintListElement(particleIterator);
				particleIterator++;
				PrintListElement(particleIterator);
				for(particleIterator = particleList.begin();particleIterator != particleList.end();++particleIterator)
				{
								PrintListElement(particleIterator);
				}

				particleIterator = particleList.begin();
				std::advance(particleIterator,1);
				std::cout << (*particleIterator)->ID << std::endl;
				
				//Particle *iterator = particleList.front();
				//std::cout << iterator -> ID << std::endl;
				//iterator = iterator -> next;

				//while(iterator != NULL)
				//{
								//std::cout << iterator -> ID << endl;
								//iterator = iterator -> next;
				//}



				return 0;
}
