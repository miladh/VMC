#ifndef ORBITALS_H
#define ORBITALS_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>

using namespace arma;
using namespace std;
using namespace libconfig;

class Orbitals
{
public:
    Orbitals();
    virtual double orbitalEvaluate(const mat &r, int qNum, int Particle) = 0;
    virtual double laplaceOrbitalEvaluate(const mat &r, int qNum, int Particle) = 0;
    virtual double getVariationalDerivative(const mat &r, int qNum, int Particle) = 0;
    virtual rowvec gradientOrbitalEvaluate(const mat &r, int qNum, int Particle) = 0;

    virtual void  setNucleusDistance()=0;
};

#endif // ORBITALS_H
