#ifndef BASICWAVEFUNCTION_H
#define BASICWAVEFUNCTION_H

#include "src/Wavefunction/wavefunction.h"

class BasicWavefunction : public Wavefunction
{
public:
    BasicWavefunction(const uint &nParticles);
    double wavefunction(const mat &r);
    double laplace(const mat &r);
    double KineticEnergy(const mat &r);
    mat gradient(const mat &r);

};

#endif // BASICWAVEFUNCTION_H
