#ifndef JASTROWWAVEFUNCTION_H
#define JASTROWWAVEFUNCTION_H

#include "src/Wavefunction/wavefunction.h"

class JastrowWavefunction : public Wavefunction
{
public:
    JastrowWavefunction();


    double jastrowFactor(int nParticles, const mat &r);
    double waveFunction(int nParticles, const mat &r);

protected:
    double rSingleParticle;
    double correlation, argument;
    double rij;
};

#endif // JASTROWWAVEFUNCTION_H
