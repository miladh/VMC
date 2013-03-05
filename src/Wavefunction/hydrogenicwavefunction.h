#ifndef HYDROGENICWAVEFUNCTION_H
#define HYDROGENICWAVEFUNCTION_H

#include "src/Wavefunction/wavefunction.h"

class HydrogenicWavefunction : public Wavefunction
{
public:
    HydrogenicWavefunction(const uint &nParticles, const int &charge);
    double wavefunction(const mat &r);
    double laplace(const mat &r);
    mat gradient(const mat &r);

private:
    int charge;
    double nFactor;

};

#endif // HYDROGENICWAVEFUNCTION_H

