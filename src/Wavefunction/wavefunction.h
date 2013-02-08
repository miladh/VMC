#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>

using namespace arma;
using namespace std;
using namespace libconfig;

class Wavefunction
{
public:
    Wavefunction();
    virtual double waveFunction(int nDimensions, int nParticles,const mat &r) = 0;
    double alpha,beta;

protected:
    double rSingleParticle;
    double correlation, argument;
    double rij;

};

#endif // WAVEFUNCTION_H
