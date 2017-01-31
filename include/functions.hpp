#ifndef _INCLUDE_FUNCTIONS_HPP_
#define _INCLUDE_FUNCTIONS_HPP_

#include "globals.hpp"
#include "classes.hpp"


void InitPositions(Particle *particle);
unsigned int NumberOfParticles();
void ComputeAccelerations(Particle* particle);

void VelocityVerlet(Particle* particle,int WRITE,FILE* output);
void VelocityVerlet(Particle* particle);

void InitVelocities(Particle* particle);
void InitVelocitiesTest(Particle* particle);
double Pressure(Particle* particle);
double* Force(Particle *particle,int i,int j);
double Distance(Particle* particle, int i, int j);

void InitBarostat(std::list<Particle*>& particles);
void InitBarostat(std::list<Particle*>& particles,std::list<Particle*>& particlesHistory);
void InitBarostat(std::vector<Particle*>& particles,std::vector<Particle*>& particleHistory);
void InitBarostat(std::vector<Particle*>& particles);
void InitBarostatNew(std::list<Particle*>& particles);
void InitBarostatFull(std::list<Particle*>& particles);
void InitBarostatFullNew(std::list<Particle*>& particles);
void InitBarostatFullNewTemp(std::list<Particle*>& particles);

void GenerateOutput(Particle *cube, std::list<Particle*> gas,FILE *combinedOut);
void GenerateOutput(Particle *cube, std::vector<Particle*> gas,FILE *combinedOut);

void CheckBoundaries(std::list<Particle*> &particles);
void CheckBoundaries(Particle* cube,std::list<Particle*> &particles);
void CheckBoundariesNew(Particle* cube,std::list<Particle*> &particles);
void CheckBoundaries(std::vector<Particle*> &particles);
void CheckBoundariesNoRemove(Particle* cube,std::list<Particle*> &particles);

void Barostat(Particle* cube, std::list<Particle*>& gas);
void Barostat(Particle* cube, std::list<Particle*>& gas,std::list<Particle*>& gasHistory);
void BarostatNew(Particle* cube, std::list<Particle*>& gas);
void BarostatNew(Particle* cube, std::list<Particle*>& gas,std::list<Particle*>& gasHistory);
void BarostatNoBoundaries(Particle* cube, std::list<Particle*>& gas);
void BarostatNoCube(std::list<Particle*>& gas);

//void BarostatTest(std::list<Particle*>& gas);
void BarostatTest(Particle* cube, std::list<Particle*>& gas);

void ComputeSoftSphere(std::list<Particle*>& gas, Particle* cube);
void ComputeSoftSphereTest(std::list<Particle*>& gas, Particle* cube);
void ComputeSoftSphereNoCube(std::list<Particle*>& gas);
void PrintAllData(Particle* cube, std::list<Particle*> gas,FILE* output);
//void eHEX(Particle* cube,FILE* output);
void eHEX(Particle* cube);
void eHEXNew(Particle* cube);
void calcPressure(Particle cube,FILE* output);
void calcTemp(Particle* cube,FILE* output);
void calcTemp(Particle* cube,FILE* output,int run);
void calcCOMTemp(double* vCOM,FILE* output);

void calcCM(Particle* particles,double *rCM, double* vCM);
void calcCM(Particle* particles,double *rCM, FILE* output);

bool fileExist(const std::string& filename);
double totErg(Particle* particles);
void writePositions(Particle* particles, std::string filename);
void writePositions(Particle* particles, std::list<Particle*>gas, std::string filename);
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
void writeHistory(std::list<Particle*>& gas, int run);
void calculateGasTemperature(std::list<Particle*> gas,FILE* output);

double calculateEnergies(Particle* cube, std::list<Particle*> gas);
double calculateEnergies(Particle* cube);
double calculateEnergiesTest(Particle* cube, std::list<Particle*> gas);
void gasStatus(std::list<Particle*> gas);

void verletBaroAccelerations(Particle* cube, std::list<Particle*> gas);
void verletBaro(Particle* cube, std::list<Particle*>& gas);
void eHEXBaro(Particle* cube, std::list<Particle*>& gas);
void eHEXBaroNewTemp(Particle* cube, std::list<Particle*>& gas);

#endif
