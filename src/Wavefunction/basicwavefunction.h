#ifndef BASICWAVEFUNCTION_H
#define BASICWAVEFUNCTION_H

#include "src/Wavefunction/wavefunction.h"

class BasicWaveFunction : public Wavefunction
{
public:
    BasicWaveFunction();
    double waveFunction(int nParticles, const mat &r);

protected:
    double rSingleParticle;
    double correlation, argument;
    double rij;
};

#endif // BASICWAVEFUNCTION_H
