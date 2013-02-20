#ifndef JASTROWWAVEFUNCTION_H
#define JASTROWWAVEFUNCTION_H

#include "src/Wavefunction/wavefunction.h"

class JastrowWavefunction : public Wavefunction
{
public:
    JastrowWavefunction();


    double jastrowFactor(int nParticles, const mat &r);
    double waveFunction(int nParticles, const mat &r);
    double laplace(int nParticles, const mat &r);

protected:
    double rSingleParticle;
    double correlation, argument;
    double rij;
private:
    int charge;
    double r1, r2, rij,r1r2;
    double E_L1,eIntEnergy,eContributor;
};

#endif // JASTROWWAVEFUNCTION_H
