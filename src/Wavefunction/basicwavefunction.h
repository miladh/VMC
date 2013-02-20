#ifndef BASICWAVEFUNCTION_H
#define BASICWAVEFUNCTION_H

#include "src/Wavefunction/wavefunction.h"

class BasicWavefunction : public Wavefunction
{
public:
    BasicWavefunction();
    double waveFunction(int nParticles, const mat &r);
    double laplace(int nParticles, const mat &r, Config *cfg);
    double KineticEnergy(int nParticles, const mat &r);

private:
    double rSingleParticle;
    double correlation, argument;
    int charge;
    double r1, r2, rij;
};

#endif // BASICWAVEFUNCTION_H
