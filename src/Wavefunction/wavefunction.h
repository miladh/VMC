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
    double waveFunction(int nDimensions, int nParticles,const mat &r, const double &alpha,const double &beta);

private:
    double rSingleParticle;
    double correlation, argument;
    double rij;


    double jastrowFactor(int nDimensions, int nParticles,const mat &r,const double &beta);
};

#endif // WAVEFUNCTION_H
