#ifndef JASTROWWAVEFUNCTION_H
#define JASTROWWAVEFUNCTION_H

#include "src/Wavefunction/wavefunction.h"

class JastrowWaveFunction : public Wavefunction
{
public:
    JastrowWaveFunction();


    double jastrowFactor(int nDimensions, int nParticles,const mat &r);
    double waveFunction(int nDimensions, int nParticles,const mat &r);
};

#endif // JASTROWWAVEFUNCTION_H
