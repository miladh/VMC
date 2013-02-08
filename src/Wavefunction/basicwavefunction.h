#ifndef BASICWAVEFUNCTION_H
#define BASICWAVEFUNCTION_H

#include "src/Wavefunction/wavefunction.h"

class BasicWaveFunction : public Wavefunction
{
public:
    BasicWaveFunction();
    double waveFunction(int nDimensions, int nParticles,const mat &r);
};

#endif // BASICWAVEFUNCTION_H
