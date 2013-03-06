#ifndef HLIKEWAVEFUNCTION_H
#define HLIKEWAVEFUNCTION_H


#include <armadillo>
#include <iostream>
#include <libconfig.h++>
#include "wavefunction.h"
#include "src/slater/slater.h"
#include "src/Jastrow/jastrow.h"

using namespace arma;
using namespace std;
using namespace libconfig;


class HLikeWavefunction : public Wavefunction
{
public:
    HLikeWavefunction(const uint &nParticles);

    double jastrowFactor(const mat &r);
    double wavefunction(const mat &r);
    double laplace(const mat &r);
    mat gradient(const mat &r);

private:
    mat dHydrogenic,dJastrow;
    uint nParticles;
};

#endif // HLIKEWAVEFUNCTION_H
