#ifndef _INCLUDE_FUNCTIONS_HPP_
#define _INCLUDE_FUNCTIONS_HPP_

#include "globals.hpp"
#include "classes.hpp"


void InitPositions(Particle *particle);
unsigned int NumberOfParticles();
void ComputeAccelerations(Particle* particle);
void VelocityVerlet(Particle* particle,int WRITE,FILE* output);
void VelocityVerlet(Particle* particle,FILE* output);
void InitVelocities(Particle* particle);
double Pressure(Particle* particle);
double* Force(Particle *particle,int i,int j);
double Distance(Particle* particle, int i, int j);
void InitBarostat(std::list<Particle*>& particles);
void InitBarostat(std::vector<Particle*>& particles);
void GenerateOutput(Particle *cube, std::list<Particle*> gas,FILE *combinedOut);
void GenerateOutput(Particle *cube, std::vector<Particle*> gas,FILE *combinedOut);
void CheckBoundaries(std::list<Particle*> &particles);
void CheckBoundaries(Particle* cube,std::list<Particle*> &particles);
void CheckBoundariesNew(Particle* cube,std::list<Particle*> &particles);
void CheckBoundaries(std::vector<Particle*> &particles);
void Barostat(Particle* cube, std::list<Particle*>& gas);
void BarostatNew(Particle* cube, std::list<Particle*>& gas);
void ComputeSoftSphere(std::list<Particle*>& gas, Particle* cube);
void PrintAllData(Particle* cube, std::list<Particle*> gas,FILE* output);
//void eHEX(Particle* cube,FILE* output);
void eHEX(Particle* cube);
void calcTemp(Particle* cube,FILE* output);
void calcCM(Particle* particles,double *rCM, double* vCM);
void calcCM(Particle* particles,double *rCM, FILE* output);
bool fileExist(const std::string& filename);
double totErg(Particle* particles);
void writePositions(Particle* particles, std::string filename);
void readPositions(Particle* particles, std::string filename);
void harmonicTrap(double* rCM, double* vCM, double* pos, Particle* particles);
void GenerateOutput(Particle *cube, std::list<Particle*> gas,int Run);
void trackParticle(Particle* cube, std::list<Particle*> gas,int partID, FILE* output);
template <typename T>
std::string numberToString(T number)
{
    std::stringstream convert;
    convert << number;
    return convert.str();
}
void printDistances(Particle* cube, std::list<Particle*> &particles,int Run);
std::string DateToString();
void rescaleVelocities(Particle* cube);
void setValues(double temp, double dq, double Eps, double Pressure, double ambienttemp);


#endif
