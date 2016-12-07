#ifndef _INCLUDE_CLASSES_H_
#define _INCLUDE_CLASSES_H_

#include "globals.hpp"

//class definitions
class Particle
{
    public:
     double *r;
     double *v;
     double *a;
     double *aOld;
     double *rnew;
     double *vnew;
     double *vhalf;
     double *vOld;
     double* r0;
     double* v0;
     double m;
     double potE;
     int type;
     int ID;
     std::string name;
     bool surface;

     //std::vector<std::pair<std::string,double> >  kineticEnergy;
     //std::vector<std::pair<std::string,double> >  potentialEnergy;
     //std::vector<std::pair<std::string,double> >  totalEnergy;

     Particle();
     Particle(std::string name);
     Particle(double x,double y, double z, double v1, double v2, double v3);
     Particle(double x,double y, double z, double v1, double v2, double v3, std::string Name);
     ~Particle();
};

#endif
